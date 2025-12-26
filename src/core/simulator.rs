//! 模拟器，负责协调地震动强度场模拟的各个步骤


use super::gmpe;
use super::b_res_sim::BResSimulator;
use super::w_res_sim::WResSimulator;
use super::eq_source::EQSource;
use super::site::Site;
use super::geo;
use super::io;
use super::grid;
use super::utilities;
use nalgebra::{DVector, DMatrix};

/// 模拟器主类
pub struct Simulator {
    eq_source: EQSource,
    /// 场地列表
    sites: Vec<Site>,
    gmpe_model_name: String,
    pub periods: Vec<f32>,
    /// 事件内残差模拟时 PCA 分量数量
    n_pcs: usize,
    /// 如果场地数量超过grid_threshold则启用网格模拟
    grid_threshold: usize,
    /// 网格模拟时的网格边长 (km)，None 表示使用默认值 0.5 km
    grid_spacing_km: Option<f32>,
    /// 输出目录
    output_dir: String,
}

impl Simulator {
    /// 创建新的模拟器实例
    pub fn new(eq_source: EQSource, sites: Vec<Site>, gmpe_model_name: String, n_pcs: usize) -> Self {
        // 定义模拟周期 (这里使用常用的周期列表，也可以根据需求调整)
        let periods: Vec<f32> = vec![
            0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0
        ];
        Self { 
            eq_source, 
            sites, 
            gmpe_model_name, 
            periods, 
            n_pcs, 
            grid_threshold: 500, 
            grid_spacing_km: None,
            output_dir: "output".to_string(),
        }
    }

    /// 设置启用网格模拟的场地数量阈值（默认 500）
    pub fn with_grid_threshold(mut self, threshold: usize) -> Self {
        self.grid_threshold = threshold;
        self
    }

    /// 设置网格间距 (km)
    pub fn with_grid_spacing(mut self, spacing: f32) -> Self {
        self.grid_spacing_km = Some(spacing);
        self
    }

    /// 设置输出目录
    pub fn with_output_dir(mut self, output_dir: String) -> Self {
        self.output_dir = output_dir;
        self
    }

    /// 运行模拟流程
    pub fn run(&self) {
        println!("\n=== 开始地震动强度模拟 ===");
        println!("震源震级: M{}", self.eq_source.m);
        println!("模拟次数: {}", self.eq_source.n_sim);
        println!("场地数量: {}", self.sites.len());
        println!("GMPE模型: {}", self.gmpe_model_name);

        // 确保输出目录存在
        if !std::path::Path::new(&self.output_dir).exists() {
            std::fs::create_dir_all(&self.output_dir).expect("Failed to create output directory");
        }

        // 0. 补全并保存震源参数
        let used_eq_source = self.eq_source.clone();
        let eq_json_path = std::path::Path::new(&self.output_dir).join("used_eq_source.json");
        let eq_json_file = std::fs::File::create(&eq_json_path).expect("Failed to create used_eq_source.json");
        serde_json::to_writer_pretty(eq_json_file, &used_eq_source).expect("Failed to write used_eq_source.json");
        println!("已保存使用的震源参数至: {:?}", eq_json_path);

        // 1. 初始化 GMPE 模型
        let gmpe_model = gmpe::create_gmpe_model(&self.gmpe_model_name)
            .expect(&format!("未知的 GMPE 模型名称: {}", self.gmpe_model_name));

        // 2. 遍历所有场地，计算中值 (Median) 和标准差 (Tau, Phi)
        println!("\n[Step 1/4] 计算 GMPE 中值和标准差...");
        let use_grid = self.sites.len() > self.grid_threshold;
        let (sim_sites, sim_coords_xy): (Vec<Site>, (Vec<f32>, Vec<f32>)) = if use_grid {
            // 生成网格点并构造网格场地（属性通过 Delaunay 插值自原始场地）
            // 使用默认 0.5 km 的网格边长；如需可配置，可在此传入 Some(value)
            let (grid_points, _nx, _ny) = grid::generate_grid_points_from_sites(&self.sites, self.grid_spacing_km);
            let sim_sites = grid::interpolate_sites_to_points(&self.sites, &grid_points);

            let mut gx = Vec::with_capacity(grid_points.len());
            let mut gy = Vec::with_capacity(grid_points.len());

            for (lon, lat) in grid_points.iter() {
                let (x, y) = geo::latlon2xy(*lon, *lat, self.eq_source.lon_0, self.eq_source.lat_0);
                gx.push(x as f32);
                gy.push(y as f32);
            }

            (sim_sites, (gx, gy))
        } else {
            // 直接使用原始场地
            let mut x_coords = Vec::with_capacity(self.sites.len());
            let mut y_coords = Vec::with_capacity(self.sites.len());
            for s in &self.sites {
                let (x, y) = geo::latlon2xy(s.lon, s.lat, self.eq_source.lon_0, self.eq_source.lat_0);
                x_coords.push(x as f32);
                y_coords.push(y as f32);
            }
            (self.sites.clone(), (x_coords, y_coords))
        };

        let n_sites = sim_sites.len();
        let n_periods = self.periods.len();
        let mut median_matrix = DMatrix::zeros(n_sites, n_periods);
        let mut phi_matrix = DMatrix::zeros(n_sites, n_periods);
        let mut tau_matrix = DMatrix::zeros(n_sites, n_periods);

        for (i, site) in sim_sites.iter().enumerate() {
            for (j, &t) in self.periods.iter().enumerate() {
                let mut site_t = site.clone();
                site_t.period1 = t as f32;
                if let Ok(res) = gmpe_model.calc(&self.eq_source, &site_t) {
                    median_matrix[(i, j)] = res.psa_median as f32;
                    phi_matrix[(i, j)] = res.psa_phi as f32;
                    tau_matrix[(i, j)] = res.psa_tau as f32;
                }
            }
        }

        // 3. 模拟事件间残差 (Between-event residuals)
        println!("\n[Step 2/4] 模拟事件间残差 (B_res)...");
        let b_res_simulator = BResSimulator::new(self.periods.clone());
        // 传入 tau_matrix 和 phi_matrix，模拟器内部会计算参考相关性并进行场地特定的缩放
        let b_res_results = b_res_simulator.simulate(&tau_matrix, &phi_matrix, self.eq_source.n_sim);

        // 4. 模拟事件内残差 (Within-event residuals)
        println!("\n[Step 3/4] 模拟事件内残差 (W_res)...");
        
        // 准备坐标数据
        let x_vec = DVector::from_vec(sim_coords_xy.0.clone());
        let y_vec = DVector::from_vec(sim_coords_xy.1.clone());
        
        // 初始化模拟器
        let mut w_res_sim = WResSimulator::new(Some(self.n_pcs), self.eq_source.seed);
        
        // 注册场地
        w_res_sim.register_sites(&x_vec, &y_vec);
        
        let periods_dvec = DVector::from_vec(self.periods.clone());
        
        // 模拟并使用预计算好的 phi_matrix 进行缩放
        let w_res_results = w_res_sim.simulate_residuals(&periods_dvec, self.eq_source.n_sim, Some(&phi_matrix));

        // 5. 组合结果并输出（包含对场地指定周期的插值与单独保存）
        println!("\n[Step 4/4] 组合结果并生成输出文件...");

        let calculated_ims = if self.eq_source.ifmedian {
            vec![median_matrix.clone(); self.eq_source.n_sim]
        } else {
            Self::combine_results(
                n_sites,
                n_periods,
                &median_matrix,
                &b_res_results,
                &w_res_results,
                self.eq_source.n_sim,
            )
        };

        let (total_ims, site_period_results) = if !use_grid {
            let total_ims = calculated_ims;
            let site_period_results = Self::interpolate_to_site_periods(&self.sites, &total_ims, &self.periods);
            (total_ims, site_period_results)
        } else {
            // 先在网格上组合
            let grid_total_ims = calculated_ims;

            // 将网格 IM 插值到原始场地上
            let site_total_ims = self.interpolate_grid_to_sites(&grid_total_ims, &sim_coords_xy, &self.sites);

            // 对站点特定周期进行插值并构造 site_period_results
            let site_period_results = Self::interpolate_to_site_periods(&self.sites, &site_total_ims, &self.periods);

            (site_total_ims, site_period_results)
        };

        // 保存所有周期的结果
        let output_dir = &self.output_dir;
        match io::save_simulation_results(output_dir, &self.periods, &total_ims) {
            Ok(_) => println!("  - 全部周期结果已保存至 {} 目录", output_dir),
            Err(e) => println!("  - 保存全部周期结果失败: {}", e),
        }

        // 保存场地文件指定周期的插值结果至单独文件
        let site_periods: Vec<f32> = self.sites.iter().map(|s| s.period1).collect();
        match io::save_site_period_results(output_dir, &site_periods, &site_period_results) {
            Ok(_) => println!("  - 场地指定周期插值结果已单独保存"),
            Err(e) => println!("  - 保存场地指定周期插值结果失败: {}", e),
        }

        println!("\n模拟完成!");
    }
}

impl Simulator {
    /// 组合中值与残差得到总 IM
    /// # 参数：
    /// - `n_sites`: 场地数量
    /// - `n_periods`: 周期数量
    /// - `median_matrix`: 中值矩阵 [n_sites x n_periods]
    /// - `b_res_results`: 事件间残差结果数组，[n_sims]，每个元素为矩阵 [n_sites x n_periods]
    /// - `w_res_results`: 事件内残差结果数组，[n_sims]，每个元素为矩阵 [n_sites x n_periods]
    /// - `n_sim`: 模拟次数
    /// # 返回：
    /// - `total_ims`: [n_sims]，每次模拟一个矩阵 [n_sites x n_periods]
    fn combine_results(
        n_sites: usize,
        n_periods: usize,
        median_matrix: &DMatrix<f32>,
        b_res_results: &[DMatrix<f32>],
        w_res_results: &[DMatrix<f32>],
        n_sim: usize,
    ) -> Vec<DMatrix<f32>> {
        let mut total_ims = Vec::with_capacity(n_sim);

        for k in 0..n_sim {
            let mut sim_matrix = DMatrix::zeros(n_sites, n_periods);
            let b_res = &b_res_results[k];
            let w_res = &w_res_results[k];

            // 组合为总 IM
            for i in 0..n_sites {
                for j in 0..n_periods {
                    let median = median_matrix[(i, j)];
                    sim_matrix[(i, j)] = median * (b_res[(i, j)] + w_res[(i, j)]).exp();
                }
            }
            total_ims.push(sim_matrix);
        }
        total_ims
    }

    /// 对每个场地在其指定周期处进行插值
    /// # 参数：
    /// - `sites`: 场地列表
    /// - `total_ims`: 总 IM 结果，[n_sims]，每个元素为矩阵 [n_sites x n_periods]
    /// - `periods`: 周期列表
    /// # 返回：
    /// - `site_period_results`: [n_sims]，每次模拟一个向量 [n_sites]，为各场地在其 `period1` 处的 IM（经插值）
    fn interpolate_to_site_periods(
        sites: &[Site],
        total_ims: &[DMatrix<f32>],
        periods: &[f32],
    ) -> Vec<Vec<f32>> {
        let n_sim = total_ims.len();
        let n_sites = sites.len();
        let n_periods = periods.len();
        let mut site_period_results: Vec<Vec<f32>> = Vec::with_capacity(n_sim);

        for k in 0..n_sim {
            let sim_matrix = &total_ims[k];
            let mut site_vec = vec![0.0; n_sites];

            for i in 0..n_sites {
                let target_t = sites[i].period1 as f32;
                // 收集当前行的值以便插值
                let row_vals: Vec<f32> = (0..n_periods)
                    .map(|jj| sim_matrix[(i, jj)])
                    .collect();
                site_vec[i] = utilities::interp_clamped_unsorted(target_t, periods, &row_vals);
            }
            site_period_results.push(site_vec);
        }
        site_period_results
    }

    /// 将网格上的 IM（每次模拟一个矩阵 [n_grid x n_periods]）插值到原始场地位置，返回 [n_sims] 的矩阵集合。
    /// 插值方法：对每个原始场地，寻找网格上 k=4 个最近邻点，采用 1/d 加权平均。
    /// 
    /// # 参数：
    /// - `grid_total_ims`: 网格上的总 IM 结果，[n_sims]，每个元素为矩阵 [n_grid x n_periods]
    /// - `grid_xy`: 网格点的 XY 坐标元组 (Vec<x>, Vec<y>)
    /// - `sites`: 原始场地列表
    /// # 返回：
    /// - `site_total_ims`: 场地位置的总 IM 结果，[n_sims]，每个元素为矩阵 [n_sites x n_periods]
    fn interpolate_grid_to_sites(
        &self,
        grid_total_ims: &[DMatrix<f32>],
        grid_xy: &(Vec<f32>, Vec<f32>),
        sites: &[Site],
    ) -> Vec<DMatrix<f32>> {
        let n_grid = grid_xy.0.len();
        let n_sites = sites.len();
        let n_periods = self.periods.len();
        let n_sims = grid_total_ims.len();

        // 预计算站点 XY
        let mut sx = Vec::with_capacity(n_sites);
        let mut sy = Vec::with_capacity(n_sites);
        for s in sites {
            let (x, y) = geo::latlon2xy(s.lon, s.lat, self.eq_source.lon_0, self.eq_source.lat_0);
            sx.push(x as f32);
            sy.push(y as f32);
        }

        let mut site_total_ims: Vec<DMatrix<f32>> = Vec::with_capacity(n_sims);

        for k in 0..n_sims {
            let mut mat = DMatrix::zeros(n_sites, n_periods);
            let grid_mat = &grid_total_ims[k];

            for i in 0..n_sites {
                // 找到 4 个最近邻网格点
                let mut idxs: Vec<usize> = (0..n_grid).collect();
                idxs.sort_by(|&a, &b| {
                    let da2 = (grid_xy.0[a] - sx[i]).powi(2) + (grid_xy.1[a] - sy[i]).powi(2);
                    let db2 = (grid_xy.0[b] - sx[i]).powi(2) + (grid_xy.1[b] - sy[i]).powi(2);
                    da2.partial_cmp(&db2).unwrap()
                });
                let k_neigh = 4.min(n_grid);
                let neigh = &idxs[..k_neigh];

                for j in 0..n_periods {
                    // IDW 权重
                    let mut wsum = 0.0;
                    let mut vsum = 0.0;
                    for &g in neigh {
                        let d = ((grid_xy.0[g] - sx[i]).powi(2) + (grid_xy.1[g] - sy[i]).powi(2)).sqrt();
                        let w = if d < 1e-9 { 1.0 } else { 1.0 / d };
                        wsum += w;
                        vsum += w * grid_mat[(g, j)];
                    }
                    mat[(i, j)] = if wsum > 0.0 { vsum / wsum } else { 0.0 };
                }
            }

            site_total_ims.push(mat);
        }

        site_total_ims
    }
}
