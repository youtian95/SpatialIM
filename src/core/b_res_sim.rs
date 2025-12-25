//! # 事件间残差模拟
//! 
//! ## 主要原理
//! 模拟事件间残差 $\delta B$。对于一次地震模拟，所有场地的 $\delta B(T)$相同；不同周期的 $\delta B(T)$ 服从多元正态分布。模拟步骤：
//! 1. 构建跨周期相关矩阵 $\rho_B(0,T_1,T_2)$：
//!     1. 自 Baker & Jayaram (2008) 获取跨周期总相关 $\rho_{\text{total}}(T_1,T_2)$；
//!     1. 自 Loth & Baker (2013) 获取事件内相关 $\rho_W(h,T_1,T_2)$，取 $h=0$ 得 $\rho_W(0,T_1,T_2)$；
//!     1. 对上面公式同一个场地的不同周期IM计算方差，得到下面的公式，可以根据下面这个公式计算未知的 $\rho_B(T_1,T_2)$：
//!         $$ \rho_{\text{total}}(0,T_1,T_2)\ \sigma(T_1)\ \sigma(T_2) = \rho_B(T_1,T_2)\ \tau(T_1)\ \tau(T_2) + \rho_W(0,T_1,T_2)\ \phi(T_1)\ \phi(T_2). $$
//! 2. 基于 $\rho_B(T_1,T_2)$ 构建相关矩阵 $\Sigma_B$，并进行 Cholesky 分解得到下三角矩阵 $\mathbf{L}$；
//! 3. 生成事件间残差向量：
//!     $$ \boldsymbol{\delta B} = \mathbf{L} \boldsymbol{u}, $$
//!     其中 $\boldsymbol{u}$ 为一次独立标准正态随机向量的实现。同一模拟中，所有场地共享同一 $\delta B(T)$ 向量。


use nalgebra::{DMatrix, DVector, Cholesky, SymmetricEigen};
use rand::prelude::*;
use rand_distr::StandardNormal;

/// 模拟地震事件间残差 (Between-event residual simulation)
pub struct BResSimulator {
    periods: Vec<f64>,
}

impl BResSimulator {
    pub fn new(periods: Vec<f64>) -> Self {
        Self {
            periods,
        }
    }

    /// Baker & Jayaram (2008) correlation model for total residuals
    /// 
    /// Computes the correlation of epsilons (or total residuals) between two periods T1 and T2.
    /// Ref: Baker JW, Jayaram N. Correlation of spectral acceleration values from NGA ground motion models. Earthquake Spectra. 2008;24:299-317.
    fn rho_total_baker_jayaram_2008(t1: f64, t2: f64) -> f64 {
        let t_min = t1.min(t2);
        let t_max = t1.max(t2);
        let pi = std::f64::consts::PI;

        let c1 = 1.0 - (pi / 2.0 - (t_max / t_min.max(0.109)).ln() * 0.366).cos();

        let c2 = if t_max < 0.2 {
            1.0 - 0.105 * (1.0 - 1.0 / (1.0 + (100.0 * t_max - 5.0).exp())) 
                * (t_max - t_min) / (t_max - 0.0099)
        } else {
            0.0 // Not used if t_max >= 0.2
        };

        let c3 = if t_max < 0.109 { c2 } else { c1 };

        let c4 = c1 + 0.5 * (c3.sqrt() - c3) * (1.0 + (pi * t_min / 0.109).cos());

        if t_max <= 0.109 {
            c2
        } else if t_min > 0.109 {
            c1
        } else if t_max < 0.2 {
            c2.min(c4)
        } else {
            c4
        }
    }

    /// Loth & Baker (2013) correlation model for epsilon
    /// 
    /// Computes the spatial cross-correlation of epsilons at multiple periods.
    /// Ref: Loth C, Baker JW. A spatial cross-correlation model of spectral accelerations at multiple periods. Earthquake Eng Struc. 2013;42:397-417.
    fn rho_epsilon_loth_baker_2013(h: f64, t1: f64, t2: f64) -> f64 {
        if t1.min(t2) < 0.01 || t1.max(t2) > 10.0 {
            // In production code, we might want to clamp or log a warning.
            // For now, we proceed, but the interpolation might be out of bounds if not handled.
            // The C++ code throws an exception. Here we'll clamp for safety in interpolation.
        }
        
        let t1_clamped = t1.max(0.01).min(10.0);
        let t2_clamped = t2.max(0.01).min(10.0);

        let t_list = [0.01, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 7.5, 10.0001];
        
        // Table II. Short range coregionalization matrix, B1
        #[rustfmt::skip]
        let b1 = [
            [0.29, 0.25, 0.23, 0.23, 0.18, 0.10, 0.06, 0.06, 0.06],
            [0.25, 0.30, 0.20, 0.16, 0.10, 0.04, 0.03, 0.04, 0.05],
            [0.23, 0.20, 0.27, 0.18, 0.10, 0.03, 0.00, 0.01, 0.02],
            [0.23, 0.16, 0.18, 0.31, 0.22, 0.14, 0.08, 0.07, 0.07],
            [0.18, 0.10, 0.10, 0.22, 0.33, 0.24, 0.16, 0.13, 0.12],
            [0.10, 0.04, 0.03, 0.14, 0.24, 0.33, 0.26, 0.21, 0.19],
            [0.06, 0.03, 0.00, 0.08, 0.16, 0.26, 0.37, 0.30, 0.26],
            [0.06, 0.04, 0.01, 0.07, 0.13, 0.21, 0.30, 0.28, 0.24],
            [0.06, 0.05, 0.02, 0.07, 0.12, 0.19, 0.26, 0.24, 0.23],
        ];

        // Table III. Long range coregionalization matrix, B2
        #[rustfmt::skip]
        let b2 = [
            [0.47, 0.40, 0.43, 0.35, 0.27, 0.15, 0.13, 0.09, 0.12],
            [0.40, 0.42, 0.37, 0.25, 0.15, 0.03, 0.04, 0.00, 0.03],
            [0.43, 0.37, 0.45, 0.36, 0.26, 0.15, 0.09, 0.05, 0.08],
            [0.35, 0.25, 0.36, 0.42, 0.37, 0.29, 0.20, 0.16, 0.16],
            [0.27, 0.15, 0.26, 0.37, 0.48, 0.41, 0.26, 0.21, 0.21],
            [0.15, 0.03, 0.15, 0.29, 0.41, 0.55, 0.37, 0.33, 0.32],
            [0.13, 0.04, 0.09, 0.20, 0.26, 0.37, 0.51, 0.49, 0.49],
            [0.09, 0.00, 0.05, 0.16, 0.21, 0.33, 0.49, 0.62, 0.60],
            [0.12, 0.03, 0.08, 0.16, 0.21, 0.32, 0.49, 0.60, 0.68],
        ];

        // Table IV. Nugget effect coregionalization matrix, B3
        #[rustfmt::skip]
        let b3 = [
            [0.24, 0.22, 0.21, 0.09, -0.02, 0.01, 0.03, 0.02, 0.01],
            [0.22, 0.28, 0.20, 0.04, -0.05, 0.00, 0.01, 0.01, -0.01],
            [0.21, 0.20, 0.28, 0.05, -0.06, 0.00, 0.04, 0.03, 0.01],
            [0.09, 0.04, 0.05, 0.26, 0.14, 0.05, 0.05, 0.05, 0.04],
            [-0.02, -0.05, -0.06, 0.14, 0.20, 0.07, 0.05, 0.05, 0.05],
            [0.01, 0.00, 0.00, 0.05, 0.07, 0.12, 0.08, 0.07, 0.06],
            [0.03, 0.01, 0.04, 0.05, 0.05, 0.08, 0.12, 0.10, 0.08],
            [0.02, 0.01, 0.03, 0.05, 0.05, 0.07, 0.10, 0.10, 0.09],
            [0.01, -0.01, 0.01, 0.04, 0.05, 0.06, 0.08, 0.09, 0.09],
        ];

        // Find intervals
        let mut index1 = 0;
        let mut index2 = 0;
        for i in 0..t_list.len() - 1 {
            if t1_clamped >= t_list[i] && t1_clamped < t_list[i+1] { index1 = i; }
            if t2_clamped >= t_list[i] && t2_clamped < t_list[i+1] { index2 = i; }
        }

        // Helper for bilinear interpolation
        let interpolate = |matrix: &[[f64; 9]; 9]| -> f64 {
            let v11 = matrix[index1][index2];
            let v12 = matrix[index1][index2 + 1];
            let v21 = matrix[index1 + 1][index2];
            let v22 = matrix[index1 + 1][index2 + 1];

            let t1_l = t_list[index1];
            let t1_r = t_list[index1 + 1];
            let t2_l = t_list[index2];
            let t2_r = t_list[index2 + 1];

            // Interpolate along T1 first
            let val1 = v11 + (v21 - v11) / (t1_r - t1_l) * (t1_clamped - t1_l);
            let val2 = v12 + (v22 - v12) / (t1_r - t1_l) * (t1_clamped - t1_l);

            // Interpolate along T2
            val1 + (val2 - val1) / (t2_r - t2_l) * (t2_clamped - t2_l)
        };

        let b1_coeff = interpolate(&b1);
        let b2_coeff = interpolate(&b2);
        let b3_coeff = interpolate(&b3);

        let mut rho = b1_coeff * (-3.0 * h / 20.0).exp() + b2_coeff * (-3.0 * h / 70.0).exp();

        if h.abs() < 1e-6 {
            rho += b3_coeff;
        }

        rho
    }

    /// # 计算 $\rho_B(T_1,T_2)$ 
    /// 公式: $$ \rho_{\text{total}}(0,T_1,T_2)\ \sigma(T_1)\ \sigma(T_2) = \rho_B(T_1,T_2)\ \tau(T_1)\ \tau(T_2) + \rho_W(0,T_1,T_2)\ \phi(T_1)\ \phi(T_2). $$
    /// 
    /// # 输入
    ///  * `t1`, `t2` - 周期 T1 和 T2
    ///  * `tau1`, `tau2` - 事件间标准差 $\tau(T_1)$ 和 $\tau(T_2)$
    /// * `phi1`, `phi2` - 事件内标准差 $\phi(T_1)$ 和 $\phi(T_2)$
    /// # 返回
    ///  * 返回计算得到的 $\rho_B(T_1,T_2)$
    fn rho_eta_combined_method(
        t1: f64, t2: f64,
        tau1: f64, tau2: f64,
        phi1: f64, phi2: f64
    ) -> f64 {
        let rho_total = Self::rho_total_baker_jayaram_2008(t1, t2);
        // For between-event correlation, we consider h=0 for the epsilon part in the derivation
        let rho_eps = Self::rho_epsilon_loth_baker_2013(0.0, t1, t2);
        
        let sigma1 = (tau1.powi(2) + phi1.powi(2)).sqrt();
        let sigma2 = (tau2.powi(2) + phi2.powi(2)).sqrt();

        let num = rho_total * sigma1 * sigma2 - rho_eps * phi1 * phi2;
        let den = tau1 * tau2;

        if den.abs() < 1e-10 { 
            if (t1 - t2).abs() < 1e-5 { 1.0 } else { 0.0 }
        } else { 
            let rho = num / den;
            // Clamp to [-1, 1] to avoid numerical issues
            rho.max(-1.0).min(1.0)
        }
    }

    /// # 执行事件间残差模拟
    /// 
    /// ## 参数
    /// * `tau_matrix` - 事件间标准差矩阵 [n_sites x n_periods]
    /// * `phi_matrix` - 事件内标准差矩阵 [n_sites x n_periods]
    /// * `n_sims` - 模拟次数
    /// 
    /// ## 返回
    /// 返回一个包含 `n_sims` 个矩阵的列表。
    /// 每个矩阵维度为 [n_sites x n_periods]，包含该次模拟中每个场地的事件间残差 (B_res)。
    /// B_res(site, T) = eta(T) * tau(site, T)
    pub fn simulate(&self, tau_matrix: &DMatrix<f64>, phi_matrix: &DMatrix<f64>, n_sims: usize) -> Vec<DMatrix<f64>> {
        let n_periods = self.periods.len();
        let n_sites = tau_matrix.nrows();

        assert_eq!(tau_matrix.ncols(), n_periods, "tau 矩阵列数必须与周期数一致");
        assert_eq!(phi_matrix.ncols(), n_periods, "phi 矩阵列数必须与周期数一致");
        assert_eq!(phi_matrix.nrows(), n_sites, "phi 矩阵行数必须与场地数一致");
        
        // 1. 计算参考 tau 和 phi 向量（取所有场地的均值）
        // 用于计算 Combined Method 的相关性矩阵
        let mut tau_ref = DVector::zeros(n_periods);
        let mut phi_ref = DVector::zeros(n_periods);
        
        for j in 0..n_periods {
            tau_ref[j] = tau_matrix.column(j).mean();
            phi_ref[j] = phi_matrix.column(j).mean();
        }

        // 2. 构建事件间残差 (eta) 的相关性矩阵
        let mut correlation_matrix = DMatrix::zeros(n_periods, n_periods);

        for i in 0..n_periods {
            for j in 0..=i {
                let t1 = self.periods[i];
                let t2 = self.periods[j];
                let tau1 = tau_ref[i];
                let phi1 = phi_ref[i];
                let tau2 = tau_ref[j];
                let phi2 = phi_ref[j];

                let rho = if i == j {
                    1.0
                } else {
                    // 使用 Combined Method 计算 eta 的相关系数
                    Self::rho_eta_combined_method(t1, t2, tau1, tau2, phi1, phi2)
                };
                
                correlation_matrix[(i, j)] = rho;
                correlation_matrix[(j, i)] = rho;
            }
        }

        // 3. 对相关性矩阵进行乔斯基 (Cholesky) 分解
        let l = match Cholesky::new(correlation_matrix.clone()) {
            Some(cholesky) => cholesky.l(),
            None => {
                // 回退方案：特征值分解
                let eigen = SymmetricEigen::new(correlation_matrix.clone());
                let mut eigenvalues = eigen.eigenvalues;
                let eigenvectors = eigen.eigenvectors;

                for i in 0..n_periods {
                    if eigenvalues[i] < 0.0 {
                        eigenvalues[i] = 0.0;
                    }
                }

                // 重构分解矩阵 L = V * sqrt(D)
                let sqrt_eigenvalues = eigenvalues.map(|v| v.sqrt());
                let d_sqrt = DMatrix::from_diagonal(&sqrt_eigenvalues);
                
                eigenvectors * d_sqrt
            }
        };
        
        let mut rng = rand::rng();
        let mut results = Vec::with_capacity(n_sims);

        // 4. 生成模拟结果
        for _ in 0..n_sims {
            // 生成独立标准正态分布随机向量 Z
            let z_data: Vec<f64> = (0..n_periods).map(|_| rng.sample(StandardNormal)).collect();
            let z = DVector::from_vec(z_data);

            // 得到标准化的事件间残差 eta (长度为 n_periods)
            // eta 对于同一次事件中的所有场地是相同的（但在缩放前）
            let eta = &l * z;

            // 构建该次模拟的残差矩阵 [n_sites x n_periods]
            // B_res[i, j] = eta[j] * tau_matrix[i, j]
            let mut b_res_matrix = DMatrix::zeros(n_sites, n_periods);
            for i in 0..n_sites {
                for j in 0..n_periods {
                    b_res_matrix[(i, j)] = eta[j] * tau_matrix[(i, j)];
                }
            }
            
            results.push(b_res_matrix);
        }

        results
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    fn pearson_corr(x: &[f64], y: &[f64]) -> f64 {
        assert_eq!(x.len(), y.len());
        let n = x.len() as f64;
        let mean_x = x.iter().sum::<f64>() / n;
        let mean_y = y.iter().sum::<f64>() / n;
        let mut num = 0.0;
        let mut den_x = 0.0;
        let mut den_y = 0.0;
        for i in 0..x.len() {
            let dx = x[i] - mean_x;
            let dy = y[i] - mean_y;
            num += dx * dy;
            den_x += dx * dx;
            den_y += dy * dy;
        }
        if den_x == 0.0 || den_y == 0.0 { return 0.0; }
        num / (den_x.sqrt() * den_y.sqrt())
    }

    #[test]
    fn corr_across_sites_same_period_is_one() {
        // Periods to simulate
        let periods = vec![0.1, 0.5, 1.0, 2.0];
        let sim = BResSimulator::new(periods.clone());

        let n_sites = 5;
        let n_periods = periods.len();

        // Tau varies by site and period (positive values)
        let mut tau_matrix = DMatrix::zeros(n_sites, n_periods);
        for i in 0..n_sites {
            for j in 0..n_periods {
                tau_matrix[(i, j)] = 0.2 + 0.05 * (i as f64) + 0.03 * (j as f64);
            }
        }

        // Phi can be constant; only used to build correlation of eta
        let mut phi_matrix = DMatrix::zeros(n_sites, n_periods);
        for i in 0..n_sites {
            for j in 0..n_periods {
                phi_matrix[(i, j)] = 0.5 + 0.02 * (j as f64);
            }
        }

        let n_sims = 500; // enough for stable correlation
        let results = sim.simulate(&tau_matrix, &phi_matrix, n_sims);

        // For every period j, compute correlations across sites over simulations
        for j in 0..n_periods {
            let mut series_by_site: Vec<Vec<f64>> = vec![vec![0.0; n_sims]; n_sites];
            for s in 0..n_sims {
                let m = &results[s];
                for i in 0..n_sites {
                    series_by_site[i][s] = m[(i, j)];
                }
            }

            // Correlation between any two sites at the same period should be ~ 1.0
            for i in 0..n_sites {
                for k in (i + 1)..n_sites {
                    let r = pearson_corr(&series_by_site[i], &series_by_site[k]);
                    assert!(r > 0.999, "period {} corr(site {}, site {}) = {}", j, i, k, r);
                }
            }
        }
    }

    #[test]
    fn corr_across_periods_is_higher_for_closer_periods() {
        // Periods
        let periods = vec![0.1, 0.5, 1.0, 2.0];
        let sim = BResSimulator::new(periods.clone());

        let n_sites = 4;
        let n_periods = periods.len();

        // Use constant tau/phi across sites per period for stability
        let mut tau_matrix = DMatrix::zeros(n_sites, n_periods);
        let mut phi_matrix = DMatrix::zeros(n_sites, n_periods);
        for j in 0..n_periods {
            let tau_j = 0.25 + 0.05 * (j as f64);
            let phi_j = 0.6 + 0.02 * (j as f64);
            for i in 0..n_sites {
                tau_matrix[(i, j)] = tau_j;
                phi_matrix[(i, j)] = phi_j;
            }
        }

        // Run simulations and compute empirical correlation across periods for a fixed site (site 0)
        let n_sims = 800;
        let results = sim.simulate(&tau_matrix, &phi_matrix, n_sims);

        // Collect series per period for site 0
        let mut series_by_period: Vec<Vec<f64>> = vec![vec![0.0; n_sims]; n_periods];
        for s in 0..n_sims {
            let m = &results[s];
            for j in 0..n_periods {
                series_by_period[j][s] = m[(0, j)];
            }
        }

        // Empirical correlation matrix
        let mut empirical = DMatrix::identity(n_periods, n_periods);
        for i in 0..n_periods {
            for j in 0..n_periods {
                if i == j { continue; }
                let r = pearson_corr(&series_by_period[i], &series_by_period[j]);
                empirical[(i, j)] = r;
            }
        }

        // Check trend: correlations should generally decrease with increasing |ΔT|.
        // We enforce two mild conditions per base period: (1) nearest >= farthest (with tol);
        // (2) Pearson correlation between distance and correlation is negative enough.
        let tol = 0.08; // allow statistical fluctuation
        for i in 0..n_periods {
            // Build arrays of distances and correlations for j != i
            let mut dist_corr: Vec<(f64, f64)> = Vec::new();
            for j in 0..n_periods {
                if i == j { continue; }
                let d = (periods[i] - periods[j]).abs();
                dist_corr.push((d, empirical[(i, j)]));
            }
            // Sort by ascending distance
            dist_corr.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            let r_nearest = dist_corr.first().unwrap().1;
            let r_farthest = dist_corr.last().unwrap().1;
            assert!(r_nearest + tol >= r_farthest,
                "period i={} nearest corr {} should >= farthest {} within tol {}",
                i, r_nearest, r_farthest, tol);

            // Compute Pearson correlation between distance and correlation
            let dists: Vec<f64> = dist_corr.iter().map(|(d, _)| *d).collect();
            let cors: Vec<f64> = dist_corr.iter().map(|(_, r)| *r).collect();
            let trend = pearson_corr(&dists, &cors);
            assert!(trend < -0.2,
                "period i={} distance-correlation trend not negative enough: {} (dists {:?}, cors {:?})",
                i, trend, dists, cors);
        }
    }
}
