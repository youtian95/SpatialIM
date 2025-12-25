use super::{GMPEModel, GMPEResult, IMType};
use crate::core::eq_source::{EQSource, Region};
use crate::core::site::Site;


/// 标准周期列表
pub const PERIODS_21: [f64; 21] = [
    0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3,
    0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,
];

/// 系数结构体，用于存储 PSA(T), PGA, PGV 的值
#[derive(Debug, Clone)]
struct Coefficient {
    /// 对应 periods_21 的 PSA 系数
    psa: Vec<f64>,
    /// PGA 系数
    pga: f64,
    /// PGV 系数
    pgv: f64,
}

impl Coefficient {
    fn new(psa: Vec<f64>, pga: f64, pgv: f64) -> Self {
        Self { psa, pga, pgv }
    }

    /// 根据 IMType 获取对应的系数值
    fn get(&self, im_type: IMType) -> f64 {
        match im_type {
            IMType::PGA => self.pga,
            IMType::PGV => self.pgv,
            IMType::PSA(t) => {
                let periods = &PERIODS_21;
                // 仅尝试精确匹配，非标准周期应在 calc_log_y 中处理
                if let Some(i) = periods.iter().position(|&p| (p - t).abs() < 1e-5) {
                    self.psa[i]
                } else {
                    panic!("Coefficient::get: Unsupported or non-standard period: {}.", t);
                }
            }
        }
    }
}

/// CB14 模型实现
pub struct CB14 {
    /// 模型参数
    c: Coefficient,
    n: Coefficient,
    c0: Coefficient,
    c1: Coefficient,
    c2: Coefficient,
    c3: Coefficient,
    c4: Coefficient,
    c5: Coefficient,
    c6: Coefficient,
    c7: Coefficient,
    c8: Coefficient,
    c9: Coefficient,
    c10: Coefficient,
    c11: Coefficient,
    c12: Coefficient,
    c13: Coefficient,
    c14: Coefficient,
    c15: Coefficient,
    c16: Coefficient,
    c17: Coefficient,
    c18: Coefficient,
    c19: Coefficient,
    c20: Coefficient,
    dc20_ji: Coefficient,
    dc20_ch: Coefficient,
    k1: Coefficient,
    k2: Coefficient,
    k3: Coefficient,
    a2: Coefficient,
    h1: Coefficient,
    h2: Coefficient,
    h3: Coefficient,
    h4: Coefficient,
    h5: Coefficient,
    h6: Coefficient,
    tau1: Coefficient,
    tau2: Coefficient,
    phi1: Coefficient,
    phi2: Coefficient,
    phi_ln_af: Coefficient,
    rho_ln_pga_ln_y: Coefficient,
}

impl CB14 {
    pub fn new() -> Self {

        let c = Coefficient::new(vec![1.88; 21], 1.88, 1.88);
        let n = Coefficient::new(vec![1.18; 21], 1.18, 1.18);
        let c0 = Coefficient::new(vec![-4.365, -4.348, -4.024, -3.479, -3.293, -3.666, -4.866, -5.411, -5.962, -6.403, -7.566, -8.379, -9.841, -11.011, -12.469, -12.969, -13.306, -14.020, -14.558, -15.509, -15.975], -4.416, -2.895);
        let c1 = Coefficient::new(vec![0.977, 0.976, 0.931, 0.887, 0.902, 0.993, 1.267, 1.366, 1.458, 1.528, 1.739, 1.872, 2.021, 2.180, 2.270, 2.271, 2.150, 2.132, 2.116, 2.223, 2.132], 0.984, 1.510);
        let c2 = Coefficient::new(vec![0.533, 0.549, 0.628, 0.674, 0.726, 0.698, 0.510, 0.447, 0.274, 0.193, -0.020, -0.121, -0.042, -0.069, 0.047, 0.149, 0.368, 0.726, 1.027, 0.169, 0.367], 0.537, 0.270);
        let c3 = Coefficient::new(vec![-1.485, -1.488, -1.494, -1.388, -1.469, -1.572, -1.669, -1.750, -1.711, -1.770, -1.594, -1.577, -1.757, -1.707, -1.621, -1.512, -1.315, -1.506, -1.721, -0.756, -0.800], -1.499, -1.299);
        let c4 = Coefficient::new(vec![-0.499, -0.501, -0.517, -0.615, -0.596, -0.536, -0.490, -0.451, -0.404, -0.321, -0.426, -0.440, -0.443, -0.527, -0.630, -0.768, -0.890, -0.885, -0.878, -1.077, -1.282], -0.496, -0.453);
        let c5 = Coefficient::new(vec![-2.773, -2.772, -2.782, -2.791, -2.745, -2.633, -2.458, -2.421, -2.392, -2.376, -2.303, -2.296, -2.232, -2.158, -2.063, -2.104, -2.051, -1.986, -2.021, -2.179, -2.244], -2.773, -2.466);
        let c6 = Coefficient::new(vec![0.248, 0.247, 0.246, 0.240, 0.227, 0.210, 0.183, 0.182, 0.189, 0.195, 0.185, 0.186, 0.186, 0.169, 0.158, 0.158, 0.148, 0.135, 0.135, 0.165, 0.180], 0.248, 0.204);
        let c7 = Coefficient::new(vec![6.753, 6.502, 6.291, 6.317, 6.861, 7.294, 8.031, 8.385, 7.534, 6.990, 7.012, 6.902, 5.522, 5.650, 5.795, 6.632, 6.759, 7.978, 8.538, 8.468, 6.564], 6.768, 5.837);
        let c8 = Coefficient::new(vec![0.0; 21], 0.0, 0.0);
        let c9 = Coefficient::new(vec![-0.214, -0.208, -0.213, -0.244, -0.266, -0.229, -0.211, -0.163, -0.150, -0.131, -0.159, -0.153, -0.090, -0.105, -0.058, -0.028, 0.0, 0.0, 0.0, 0.0, 0.0], -0.212, -0.168);
        let c10 = Coefficient::new(vec![0.720, 0.730, 0.759, 0.826, 0.815, 0.831, 0.749, 0.764, 0.716, 0.737, 0.738, 0.718, 0.795, 0.556, 0.480, 0.401, 0.206, 0.105, 0.0, 0.0, 0.0], 0.720, 0.305);
        let c11 = Coefficient::new(vec![1.094, 1.149, 1.290, 1.449, 1.535, 1.615, 1.877, 2.069, 2.205, 2.306, 2.398, 2.355, 1.995, 1.447, 0.330, -0.514, -0.848, -0.793, -0.748, -0.664, -0.576], 1.090, 1.713);
        let c12 = Coefficient::new(vec![2.191, 2.189, 2.164, 2.138, 2.446, 2.969, 3.544, 3.707, 3.343, 3.334, 3.544, 3.016, 2.616, 2.470, 2.108, 1.327, 0.601, 0.568, 0.356, 0.075, -0.027], 2.186, 2.602);
        let c13 = Coefficient::new(vec![1.416, 1.453, 1.476, 1.549, 1.772, 1.916, 2.161, 2.465, 2.766, 3.011, 3.203, 3.333, 3.054, 2.562, 1.453, 0.657, 0.367, 0.306, 0.268, 0.374, 0.297], 1.420, 2.457);
        let c14 = Coefficient::new(vec![-0.0070, -0.0167, -0.0422, -0.0663, -0.0794, -0.0294, 0.0642, 0.0968, 0.1441, 0.1597, 0.1410, 0.1474, 0.1764, 0.2593, 0.2881, 0.3112, 0.3478, 0.3747, 0.3382, 0.3754, 0.3506], -0.0064, 0.1060);
        let c15 = Coefficient::new(vec![-0.207, -0.199, -0.202, -0.339, -0.404, -0.416, -0.407, -0.311, -0.172, -0.084, 0.085, 0.233, 0.411, 0.479, 0.566, 0.562, 0.534, 0.522, 0.477, 0.321, 0.174], -0.202, 0.332);
        
        let c16 = Coefficient::new(vec![0.390, 0.387, 0.378, 0.295, 0.322, 0.384, 0.417, 0.404, 0.466, 0.528, 0.540, 0.638, 0.776, 0.771, 0.748, 0.763, 0.686, 0.691, 0.670, 0.757, 0.621], 0.393, 0.585);
        let c17 = Coefficient::new(vec![0.0981, 0.1009, 0.1095, 0.1226, 0.1165, 0.0998, 0.0760, 0.0571, 0.0437, 0.0323, 0.0209, 0.0092, -0.0082, -0.0131, -0.0187, -0.0258, -0.0311, -0.0413, -0.0281, -0.0205, 0.0009], 0.0977, 0.0517);
        let c18 = Coefficient::new(vec![0.0334, 0.0327, 0.0331, 0.0270, 0.0288, 0.0325, 0.0388, 0.0437, 0.0463, 0.0508, 0.0432, 0.0405, 0.0420, 0.0426, 0.0380, 0.0252, 0.0236, 0.0102, 0.0034, 0.0050, 0.0099], 0.0333, 0.0327);
        let c19 = Coefficient::new(vec![0.00755, 0.00759, 0.00790, 0.00803, 0.00811, 0.00744, 0.00716, 0.00688, 0.00556, 0.00458, 0.00401, 0.00388, 0.00420, 0.00409, 0.00424, 0.00448, 0.00345, 0.00603, 0.00805, 0.00280, 0.00458], 0.00757, 0.00613);
        let c20 = Coefficient::new(vec![-0.0055, -0.0055, -0.0057, -0.0063, -0.0070, -0.0073, -0.0069, -0.0060, -0.0055, -0.0049, -0.0037, -0.0027, -0.0016, -0.0006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], -0.0055, -0.0017);
        
        let dc20_ji = Coefficient::new(vec![-0.0035, -0.0035, -0.0034, -0.0037, -0.0037, -0.0034, -0.0030, -0.0031, -0.0033, -0.0035, -0.0034, -0.0034, -0.0032, -0.0030, -0.0019, -0.0005, 0.0, 0.0, 0.0, 0.0, 0.0], -0.0035, -0.0006);
        let dc20_ch = Coefficient::new(vec![0.0036, 0.0036, 0.0037, 0.0040, 0.0039, 0.0042, 0.0042, 0.0041, 0.0036, 0.0031, 0.0028, 0.0025, 0.0016, 0.0006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.0036, 0.0017);
        
        let k1 = Coefficient::new(vec![865.0, 865.0, 908.0, 1054.0, 1086.0, 1032.0, 878.0, 748.0, 654.0, 587.0, 503.0, 457.0, 410.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0], 865.0, 400.0);
        let k2 = Coefficient::new(vec![-1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0.0, 0.0, 0.0, 0.0, 0.0], -1.186, -1.955);
        let k3 = Coefficient::new(vec![1.839, 1.840, 1.841, 1.843, 1.845, 1.847, 1.852, 1.856, 1.861, 1.865, 1.874, 1.883, 1.906, 1.929, 1.974, 2.019, 2.110, 2.200, 2.291, 2.517, 2.744], 1.839, 1.929);
        
        let a2 = Coefficient::new(vec![0.168, 0.166, 0.167, 0.173, 0.198, 0.174, 0.198, 0.204, 0.185, 0.164, 0.160, 0.184, 0.216, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596], 0.167, 0.596);
        
        let h1 = Coefficient::new(vec![0.242, 0.244, 0.246, 0.251, 0.260, 0.259, 0.254, 0.237, 0.206, 0.210, 0.226, 0.217, 0.154, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117], 0.241, 0.117);
        let h2 = Coefficient::new(vec![1.471, 1.467, 1.467, 1.449, 1.435, 1.449, 1.461, 1.484, 1.581, 1.586, 1.544, 1.554, 1.626, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616], 1.474, 1.616);
        let h3 = Coefficient::new(vec![-0.714, -0.711, -0.713, -0.701, -0.695, -0.708, -0.715, -0.721, -0.787, -0.795, -0.770, -0.770, -0.780, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733], -0.715, -0.733);
        let h4 = Coefficient::new(vec![1.0; 21], 0.0, 0.0);
        let h5 = Coefficient::new(vec![-0.336, -0.339, -0.338, -0.338, -0.347, -0.391, -0.449, -0.393, -0.339, -0.447, -0.525, -0.407, -0.371, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128], -0.337, -0.128);
        let h6 = Coefficient::new(vec![-0.270, -0.263, -0.259, -0.263, -0.219, -0.201, -0.099, -0.198, -0.210, -0.121, -0.086, -0.281, -0.285, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756], -0.270, -0.756);

        let tau1 = Coefficient::new(vec![0.404, 0.417, 0.446, 0.508, 0.504, 0.445, 0.382, 0.339, 0.340, 0.340, 0.356, 0.379, 0.430, 0.470, 0.497, 0.499, 0.500, 0.543, 0.534, 0.523, 0.466], 0.409, 0.317);
        let tau2 = Coefficient::new(vec![0.325, 0.326, 0.344, 0.377, 0.418, 0.426, 0.387, 0.338, 0.316, 0.300, 0.264, 0.263, 0.326, 0.353, 0.399, 0.400, 0.417, 0.393, 0.421, 0.438, 0.438], 0.322, 0.297);
        let phi1 = Coefficient::new(vec![0.734, 0.738, 0.747, 0.777, 0.782, 0.769, 0.769, 0.761, 0.744, 0.727, 0.690, 0.663, 0.606, 0.579, 0.541, 0.529, 0.527, 0.521, 0.502, 0.457, 0.441], 0.734, 0.655);
        let phi2 = Coefficient::new(vec![0.492, 0.496, 0.503, 0.520, 0.535, 0.543, 0.543, 0.552, 0.545, 0.568, 0.593, 0.611, 0.633, 0.628, 0.603, 0.588, 0.578, 0.559, 0.551, 0.546, 0.543], 0.492, 0.494);
        let phi_ln_af = Coefficient::new(vec![0.300; 21], 0.300, 0.300);
        let rho_ln_pga_ln_y = Coefficient::new(vec![1.000, 0.998, 0.986, 0.938, 0.887, 0.870, 0.876, 0.870, 0.850, 0.819, 0.743, 0.684, 0.562, 0.467, 0.364, 0.298, 0.234, 0.202, 0.184, 0.176, 0.154], 1.000, 0.684);

        CB14 { 
            c, n, 
            c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20,
            dc20_ji, dc20_ch, k1, k2, k3, a2, h1, h2, h3, h4, h5, h6,
            tau1, tau2, phi1, phi2, phi_ln_af, rho_ln_pga_ln_y
        }
    }

    /// 计算给定周期 T 下的对数地震动强度，返回 ln(Y)，可能是 PGA, PGV, PSA(T)
    fn calc_log_y(&self, eq: &EQSource, site: &Site, t: f64, im_type: IMType) -> f64 {
        // Interpolation for non-standard periods
        if let IMType::PSA(period) = im_type {
            let periods = &PERIODS_21;
            
            // Handle out of bounds (Clamp to nearest standard period)
            if period < periods[0] {
                 return self.calc_log_y(eq, site, periods[0], IMType::PSA(periods[0]));
            }
            if period > periods[periods.len() - 1] {
                 return self.calc_log_y(eq, site, periods[periods.len() - 1], IMType::PSA(periods[periods.len() - 1]));
            }

            // Check if period is standard (exact match)
            let is_standard = periods.iter().any(|&p| (p - period).abs() < 1e-5);
            
            if !is_standard {
                // Find neighbors
                let mut idx_low = 0;
                let mut idx_high = 0;
                let mut found = false;
                
                for (i, &p) in periods.iter().enumerate() {
                    if p > period {
                        idx_high = i;
                        if i > 0 {
                            idx_low = i - 1;
                            found = true;
                        }
                        break;
                    }
                }
                
                if found {
                    let t_low = periods[idx_low];
                    let t_high = periods[idx_high];
                    let log_y_low = self.calc_log_y(eq, site, t_low, IMType::PSA(t_low));
                    let log_y_high = self.calc_log_y(eq, site, t_high, IMType::PSA(t_high));
                    
                    return crate::core::utilities::linear_interp(t_low.ln(), log_y_low, t_high.ln(), log_y_high, period.ln());
                }
            }
        }

        let mag = eq.m;
        let r_rup = site.r_rup.unwrap_or_else(|| eq.calc_rrup(site));
        let r_jb = site.r_jb.unwrap_or_else(|| eq.calc_rjb(site));
        let r_x = site.r_x.unwrap_or_else(|| eq.calc_rx(site));
        let w = eq.w.unwrap();
        let dip = eq.delta();
        let z_tor = eq.z_tor.unwrap();
        let f_rv = eq.get_f_rv();
        let f_nm = eq.get_f_nm();
        let vs30 = site.vs30;
        let site_rock = Site {
            vs30: 1100.0,
            ..site.clone()
        };
        let a1100: Option<f64> = if vs30 >= 1100.0 {
            // 此时是基岩，不需要迭代计算基岩PGA a1100
            None
        } else {
            Some(self.calc_log_y(eq, &site_rock, 0.0, IMType::PGA).exp())
        };
        let sj = eq.get_if_is_japan_site();
        let zhyp = eq.zhyp.unwrap();
        let z25 = site.get_z25(sj);

        let log_y = self.f_mag(mag, im_type) 
            + self.f_dis(mag, r_rup, im_type) 
            + self.f_flt(mag, f_rv, f_nm, im_type) 
            + self.f_hng(mag, r_rup, r_jb, r_x, w, dip, z_tor, im_type) 
            + self.f_site(vs30, a1100, im_type, sj) 
            + self.f_sed(z25, sj, im_type) 
            + self.f_hyp(zhyp, mag, im_type) 
            + self.f_dip(mag, dip, im_type) 
            + self.f_atn(r_rup, eq.region, im_type);

        if matches!(im_type, IMType::PSA(_)) && t < 0.25 {
            let log_pga = self.calc_log_y(eq, site, 0.0, IMType::PGA);
            if log_y < log_pga {
                return log_pga;
            }
        }
        
        log_y
    }

    /// Magnitude Term
    fn f_mag(&self, mag: f64, im_type: IMType) -> f64 {
        let c0 = self.c0.get(im_type);
        let c1 = self.c1.get(im_type);
        let c2 = self.c2.get(im_type);
        let c3 = self.c3.get(im_type);
        let c4 = self.c4.get(im_type);
        if mag <= 4.5 {
            c0 + c1 * mag
        } else if mag <= 5.5 {
            c0 + c1 * mag + (mag - 4.5) * c2
        } else if mag <= 6.5 {
            c0 + c1 * mag + (mag - 4.5) * c2 + (mag - 5.5) * c3
        } else {
            c0 + c1 * mag + (mag - 4.5) * c2 + (mag - 5.5) * c3 + (mag - 6.5) * c4
        }
    }

    /// Geometric Attenuation Term
    fn f_dis(&self, mag: f64, r_rup: f64,  im_type: IMType) -> f64 {
        let c5 = self.c5.get(im_type);
        let c6 = self.c6.get(im_type);
        let c7 = self.c7.get(im_type);

        (c5 + c6*mag) * (r_rup.powi(2) + c7.powi(2)).sqrt().ln()
    }

    /// Faulting Style Term
    fn f_flt(&self, mag: f64, f_rv: bool, f_nm: bool, im_type: IMType) -> f64 {
        let c8 = self.c8.get(im_type);
        let c9 = self.c9.get(im_type);

        let f_rv = if f_rv { 1.0 } else { 0.0 };
        let f_nm = if f_nm { 1.0 } else { 0.0 };

        let f_flt_f = c8 * f_rv + c9 * f_nm;
        let f_flt_m = if mag <= 4.5 {
            0.0
        } else if mag <= 5.5 {
            mag - 4.5
        } else {
            1.0
        };

        f_flt_f * f_flt_m
        
    }

    /// Hanging Wall Term
    fn f_hng(&self, mag: f64, r_rup: f64, r_jb: f64, r_x: f64, w: f64, dip: f64, z_tor: f64, im_type: IMType) -> f64 {
        let c10 = self.c10.get(im_type);
        if c10 == 0.0 { return 0.0; }

        let a2 = self.a2.get(im_type);
        let h1 = self.h1.get(im_type);
        let h2 = self.h2.get(im_type);
        let h3 = self.h3.get(im_type);
        let h4 = self.h4.get(im_type);
        let h5 = self.h5.get(im_type);
        let h6 = self.h6.get(im_type);

        // f_hng_Rx
        let r1 = w * dip.to_radians().cos();
        let r2 = 62.0 * mag - 350.0;
        
        let f_hng_rx = if r_x < 0.0 {
            0.0
        } else if r_x < r1 {
            h1 + h2 * (r_x / r1) + h3 * (r_x / r1).powi(2)
        } else {
            let val = h4 + h5 * (r_x - r1) / (r2 - r1) + h6 * ((r_x - r1) / (r2 - r1)).powi(2);
            val.max(0.0)
        };

        // f_hng_Rrup
        let f_hng_rrup = if r_rup == 0.0 {
            1.0
        } else {
            (r_rup - r_jb) / r_rup
        };

        // f_hng_M
        let f_hng_m = if mag <= 5.5 {
            0.0
        } else if mag <= 6.5 {
            (mag - 5.5) * (1.0 + a2 * (mag - 6.5))
        } else {
            1.0 + a2 * (mag - 6.5)
        };

        // f_hng_Z
        let f_hng_z = if z_tor <= 16.66 {
            1.0 - 0.06 * z_tor
        } else {
            0.0
        };

        // f_hng_delta
        let f_hng_delta = (90.0 - dip) / 45.0;

        c10 * f_hng_rx * f_hng_rrup * f_hng_m * f_hng_z * f_hng_delta
    }

    /// Site Term
    /// Parameters:
    /// - vs30: Shear-wave velocity in the top 30 meters (m/s)
    /// - a1100: PGA on rock (vs30 = 1100 m/s) in g
    /// - im_type: Intensity Measure Type (PGA, PGV, PSA(T))
    /// - sj: Japan site condition flag. If true, apply Japan-specific site term adjustments.
    fn f_site(&self, vs30: f64, a1100: Option<f64>, im_type: IMType, sj: bool) -> f64 {
        let c11 = self.c11.get(im_type);
        let k1 = self.k1.get(im_type);
        let k2 = self.k2.get(im_type);
        let n = self.n.get(im_type);
        let c = self.c.get(im_type);

        // f_site_G
        let f_site_g = if vs30 <= k1 {
            let term1 = c11 * (vs30 / k1).ln();
            let term2 = k2 * ((a1100.unwrap() + c * (vs30 / k1).powf(n)).ln() - (a1100.unwrap() + c).ln());
            term1 + term2
        } else {
            (c11 + k2 * n) * (vs30 / k1).ln()
        };

        // f_site_J (Japan)
        let f_site_j = if sj {
            let c12 = self.c12.get(im_type);
            let c13 = self.c13.get(im_type);
            
            if vs30 <= 200.0 {
                (c12 + k2 * n) * ((vs30 / k1).ln() - (200.0 / k1).ln())
            } else {
                (c13 + k2 * n) * (vs30 / k1).ln()
            }
        } else {
            0.0
        };

        f_site_g + f_site_j
    }

    /// Sediment Depth Term
    fn f_sed(&self, z25: f64, sj: bool, im_type: IMType) -> f64 {
        let c14 = self.c14.get(im_type);
        let c15 = self.c15.get(im_type);
        let c16 = self.c16.get(im_type);
        let k3 = self.k3.get(im_type);

        let sj_val = if sj { 1.0 } else { 0.0 };

        if z25 <= 1.0 {
            (c14 + c15 * sj_val) * (z25 - 1.0)
        } else if z25 <= 3.0 {
            0.0
        } else {
            c16 * k3 * (-0.75f64).exp() * (1.0 - (-0.25 * (z25 - 3.0)).exp())
        }
    }

    /// Hypocenter Depth Term
    fn f_hyp(&self, zhyp: f64, mag: f64, im_type: IMType) -> f64 {
        let c17 = self.c17.get(im_type);
        let c18 = self.c18.get(im_type);

        let f_hyp_h = if zhyp <= 7.0 {
            0.0
        } else if zhyp <= 20.0 {
            zhyp - 7.0
        } else {
            13.0
        };

        let f_hyp_m = if mag <= 5.5 {
            c17
        } else if mag <= 6.5 {
            c17 + (c18 - c17) * (mag - 5.5)
        } else {
            c18
        };

        f_hyp_h * f_hyp_m
    }

    /// Fault Dip Term
    fn f_dip(&self, mag: f64, dip: f64, im_type: IMType) -> f64 {
        let c19 = self.c19.get(im_type);
        
        if mag <= 4.5 {
            c19 * dip
        } else if mag <= 5.5 {
            c19 * (5.5 - mag) * dip
        } else {
            0.0
        }
    }

    /// Anelastic Attenuation Term
    fn f_atn(&self, r_rup: f64, region: Region, im_type: IMType) -> f64 {
        if r_rup <= 80.0 {
            return 0.0;
        }

        let c20 = self.c20.get(im_type);
        let dc20 = match region {
            Region::ChinaTurkey => self.dc20_ch.get(im_type), // China/Turkey
            Region::Italy | Region::Japan => self.dc20_ji.get(im_type), // Italy or Japan
            _ => 0.0,
        };

        (c20 + dc20) * (r_rup - 80.0)
    }
    
    /// standard deviations
    pub fn calc_std_log_y(&self, eq: &EQSource, site: &Site, im_type: IMType) -> (f64, f64, f64) {
        // Interpolation for non-standard periods
        if let IMType::PSA(period) = im_type {
            let periods = &PERIODS_21;
            
            // Handle out of bounds (Clamp to nearest standard period)
            if period < periods[0] {
                 return self.calc_std_log_y(eq, site, IMType::PSA(periods[0]));
            }
            if period > periods[periods.len() - 1] {
                 return self.calc_std_log_y(eq, site, IMType::PSA(periods[periods.len() - 1]));
            }

            // Check if period is standard (exact match)
            let is_standard = periods.iter().any(|&p| (p - period).abs() < 1e-5);
            
            if !is_standard {
                // Find neighbors
                let mut idx_low = 0;
                let mut idx_high = 0;
                let mut found = false;
                
                for (i, &p) in periods.iter().enumerate() {
                    if p > period {
                        idx_high = i;
                        if i > 0 {
                            idx_low = i - 1;
                            found = true;
                        }
                        break;
                    }
                }
                
                if found {
                    let t_low = periods[idx_low];
                    let t_high = periods[idx_high];
                    let (sigma_low, tau_low, phi_low) = self.calc_std_log_y(eq, site, IMType::PSA(t_low));
                    let (sigma_high, tau_high, phi_high) = self.calc_std_log_y(eq, site, IMType::PSA(t_high));
                    
                    let sigma = crate::core::utilities::linear_interp(t_low.ln(), sigma_low, t_high.ln(), sigma_high, period.ln());
                    let tau = crate::core::utilities::linear_interp(t_low.ln(), tau_low, t_high.ln(), tau_high, period.ln());
                    let phi = crate::core::utilities::linear_interp(t_low.ln(), phi_low, t_high.ln(), phi_high, period.ln());
                    return (sigma, tau, phi);
                }
            }
        }

        let mag = eq.m;
        let vs30 = site.vs30;
        
        // Calculate A1100
        let site_rock = Site { vs30: 1100.0, ..site.clone() };
        let ln_a1100 = self.calc_log_y(eq, &site_rock, 0.0, IMType::PGA);
        let a1100 = ln_a1100.exp();

        // Get Tau_lnY and Phi_lnY for current IM
        let tau_ln_y = self.get_tau_ln_y(mag, im_type);
        let phi_ln_y = self.get_phi_ln_y(mag, im_type);
        let phi_ln_af = self.phi_ln_af.get(im_type);
        
        // Phi_lnYb
        // phi_ln_y_b = sqrt(phi_ln_y^2 - phi_ln_af^2)
        // Ensure non-negative
        let phi_ln_y_b = (phi_ln_y.powi(2) - phi_ln_af.powi(2)).max(0.0).sqrt();
        let tau_ln_y_b = tau_ln_y;

        // Get Tau_lnPGA and Phi_lnPGA
        let (tau_ln_pga_b, phi_ln_pga_b) = if im_type == IMType::PGA {
            (tau_ln_y_b, phi_ln_y_b)
        } else {
            let tau_ln_pga = self.get_tau_ln_y(mag, IMType::PGA);
            let phi_ln_pga = self.get_phi_ln_y(mag, IMType::PGA);
            let phi_ln_af_pga = self.phi_ln_af.get(IMType::PGA);
            let phi_ln_pga_b = (phi_ln_pga.powi(2) - phi_ln_af_pga.powi(2)).max(0.0).sqrt();
            (tau_ln_pga, phi_ln_pga_b)
        };

        let rho = self.rho_ln_pga_ln_y.get(im_type);
        let alpha = self.calc_alpha(vs30, a1100, im_type);

        // Eq 29
        let tau_sq = tau_ln_y_b.powi(2) + alpha.powi(2) * tau_ln_pga_b.powi(2) + 2.0 * alpha * rho * tau_ln_y_b * tau_ln_pga_b;
        let tau = tau_sq.sqrt();

        // Eq 30
        let phi_sq = phi_ln_y_b.powi(2) + phi_ln_af.powi(2) + alpha.powi(2) * phi_ln_pga_b.powi(2) + 2.0 * alpha * rho * phi_ln_y_b * phi_ln_pga_b;
        let phi = phi_sq.sqrt();

        let sigma = (tau.powi(2) + phi.powi(2)).sqrt();

        (sigma, tau, phi)
    }

    fn get_tau_ln_y(&self, mag: f64, im_type: IMType) -> f64 {
        let tau1 = self.tau1.get(im_type);
        let tau2 = self.tau2.get(im_type);
        if mag <= 4.5 {
            tau1
        } else if mag <= 5.5 {
            tau2 + (tau1 - tau2) * (5.5 - mag)
        } else {
            tau2
        }
    }

    fn get_phi_ln_y(&self, mag: f64, im_type: IMType) -> f64 {
        let phi1 = self.phi1.get(im_type);
        let phi2 = self.phi2.get(im_type);
        if mag <= 4.5 {
            phi1
        } else if mag <= 5.5 {
            phi2 + (phi1 - phi2) * (5.5 - mag)
        } else {
            phi2
        }
    }

    fn calc_alpha(&self, vs30: f64, a1100: f64, im_type: IMType) -> f64 {
        let k1 = self.k1.get(im_type);
        if vs30 >= k1 {
            return 0.0;
        }
        let k2 = self.k2.get(im_type);
        let c = self.c.get(im_type);
        let n = self.n.get(im_type);
        
        let term1 = 1.0 / (a1100 + c * (vs30 / k1).powf(n));
        let term2 = 1.0 / (a1100 + c);
        
        k2 * a1100 * (term1 - term2)
    }
}

impl GMPEModel for CB14 {
    fn calc(&self, eq: &EQSource, site: &Site) -> Result<GMPEResult, String> {
        let t = site.period1;
        
        // 1. Calculate PGA
        let ln_pga = self.calc_log_y(eq, site, 0.0, IMType::PGA);
        let pga_median = ln_pga.exp();
        let (pga_sigma, pga_tau, pga_phi) = self.calc_std_log_y(eq, site, IMType::PGA);

        // 2. Calculate PGV
        let ln_pgv = self.calc_log_y(eq, site, -1.0, IMType::PGV);
        let pgv_median = ln_pgv.exp();
        let (pgv_sigma, pgv_tau, pgv_phi) = self.calc_std_log_y(eq, site, IMType::PGV);

        // 3. Calculate PSA (at site period)
        let (psa_median, psa_sigma, psa_tau, psa_phi) = if t > 0.0 {
            let val = self.calc_log_y(eq, site, t, IMType::PSA(t));
            let (sigma, tau, phi) = self.calc_std_log_y(eq, site, IMType::PSA(t));
            (val.exp(), sigma, tau, phi)
        } else if t == 0.0 {
            (pga_median, pga_sigma, pga_tau, pga_phi)
        } else {
            // For PGV or invalid period, PSA is not really defined.
            (0.0, 0.0, 0.0, 0.0)
        };

        Ok(GMPEResult {
            psa_median,
            psa_sigma,
            psa_tau,
            psa_phi,
            pga_median,
            pga_sigma,
            pga_tau,
            pga_phi,
            pgv_median,
            pgv_sigma,
            pgv_tau,
            pgv_phi,
        })
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    /// 需要自己对比查看 CB14 论文中的 Figure 6 数据
    #[test]
    fn test_generate_figure_6_data() {
        use std::fs::File;
        use std::io::Write;

        let cb14 = CB14::new();
        
        // Table C1 Data
        struct Scenario {
            m: f64,
            z_bot: f64,
            w: f64,
            z_tor: f64,
            z_hyp: f64,
        }

        let scenarios = vec![
            Scenario { m: 3.5, z_bot: 15.0, w: 0.5119, z_tor: 7.1449, z_hyp: 7.5626 },
            Scenario { m: 4.5, z_bot: 15.0, w: 1.6572, z_tor: 7.1449, z_hyp: 8.2623 },
            Scenario { m: 5.5, z_bot: 15.0, w: 5.3653, z_tor: 4.2887, z_hyp: 7.2779 },
            Scenario { m: 6.5, z_bot: 15.0, w: 14.1259, z_tor: 0.8741, z_hyp: 8.8705 },
            Scenario { m: 7.5, z_bot: 15.0, w: 15.0000, z_tor: 0.0, z_hyp: 10.2267 },
        ];

        // Distances to calculate (R_rup)
        // Log-spaced roughly from 0.1 to 300
        let mut r_rups = vec![];
        let mut r = 0.1f64;
        while r <= 800.0 {
            r_rups.push(r);
            r *= 1.1; // Step factor
        }

        let mut path_pga = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path_pga.push("tests/test_result_cb14_pga_vs_rrup_StrikeSlip.csv");
        let mut file_pga = File::create(path_pga).unwrap();

        let mut path_psa = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path_psa.push("tests/test_result_cb14_psa1_vs_rrup_StrikeSlip.csv");
        let mut file_psa = File::create(path_psa).unwrap();

        // Write Headers
        writeln!(file_pga, "R_rup,M3.5,M4.5,M5.5,M6.5,M7.5").unwrap();
        writeln!(file_psa, "R_rup,M3.5,M4.5,M5.5,M6.5,M7.5").unwrap();

        for r_rup in r_rups {
            let mut row_pga = format!("{:.4}", r_rup);
            let mut row_psa = format!("{:.4}", r_rup);

            for s in &scenarios {
                // Check physical limit (closest distance is Z_TOR)
                if r_rup < s.z_tor {
                    row_pga.push_str(",");
                    row_psa.push_str(",");
                    continue;
                }

                // For vertical strike-slip:
                // R_rup^2 = R_jb^2 + Z_tor^2
                // R_jb = sqrt(R_rup^2 - Z_tor^2)
                let r_jb = (r_rup.powi(2) - s.z_tor.powi(2)).sqrt();
                let r_x = r_jb; // For vertical fault, Rx is horizontal distance

                let eq = EQSource::new(
                    true,
                    s.m,
                    1,
                    0,
                    0.0, 0.0, // Epicenter
                    Some(s.w),
                    None, // Length unknown
                    (1.0, 0.0, 0.0), // Normal along X (East) -> Strike North (Y), Dip 90.
                    0.0, // Strike-slip
                    false, // No HW
                    Some(s.z_hyp),
                    Some(s.z_tor),
                    Some(s.z_bot),
                    Region::California, // Default for NGA
                );

                let site = Site::new(
                    0,
                    0.0, 0.0, // Dummy coordinates
                    0.0,
                    1.0, // Period 1.0s
                    760.0,
                    Some(0.6068),
                    false // Not Japan
                ).with_distances(Some(r_rup), Some(r_jb), Some(r_x));

                let res_pga = cb14.calc_log_y(&eq, &site, 0.0, IMType::PGA).exp();
                let res_psa = cb14.calc_log_y(&eq, &site, 1.0, IMType::PSA(1.0)).exp();

                row_pga.push_str(&format!(",{:.6}", res_pga));
                row_psa.push_str(&format!(",{:.6}", res_psa));
            }
            
            writeln!(file_pga, "{}", row_pga).unwrap();
            writeln!(file_psa, "{}", row_psa).unwrap();
        }
    }

    #[test]
    fn test_generate_figure_7_data() {
        use std::fs::File;
        use std::io::Write;

        let cb14 = CB14::new();
        
        // Table C2 Data
        struct Scenario {
            m: f64,
            z_bot: f64,
            w: f64,
            z_tor: f64,
            z_hyp: f64,
        }

        let scenarios = vec![
            Scenario { m: 3.5, z_bot: 15.0, w: 0.5119, z_tor: 7.3116, z_hyp: 7.6374 },
            Scenario { m: 4.5, z_bot: 15.0, w: 1.6572, z_tor: 7.3116, z_hyp: 8.3663 },
            Scenario { m: 5.5, z_bot: 15.0, w: 5.3653, z_tor: 7.3116, z_hyp: 10.3008 },
            Scenario { m: 6.5, z_bot: 15.0, w: 16.0763, z_tor: 3.6324, z_hyp: 11.6288 },
            Scenario { m: 7.5, z_bot: 15.0, w: 20.5595, z_tor: 0.4622, z_hyp: 10.6889 },
        ];

        let mut path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("tests/test_result_cb14_pga_vs_rrup_ReverseFault.csv");
        let mut file = File::create(path).unwrap();
        
        // Header
        writeln!(file, "R_rup,M3.5,M4.5,M5.5,M6.5,M7.5").unwrap();

        // Log-spaced R_rup from 0.1 to 800
        let mut r_rups = vec![];
        let mut r = 0.1f64;
        while r <= 800.0 {
            r_rups.push(r);
            r *= 1.1; 
        }

        for r_rup in r_rups {
            let mut row = format!("{:.4}", r_rup);
            
            for s in &scenarios {
                // Check physical limit (closest distance is Z_TOR)
                if r_rup < s.z_tor {
                    row.push_str(",");
                    continue;
                }

                // Analytical calculation of Rx from R_rup for Reverse Fault (45 deg dip, HW)
                // Geometry:
                // Region 1: Site is "above" the top edge extension. Closest point is top edge.
                //           R_rup^2 = Rx^2 + Z_tor^2
                // Region 2: Site is "above" the fault plane. Closest point is on the plane.
                //           R_rup = Rx * sin(delta) + Z_tor * cos(delta)
                // Region 3: Site is "beyond" the bottom edge extension. Closest point is bottom edge.
                //           R_rup^2 = (Rx - W*cos(delta))^2 + Z_rup_bot^2
                
                let cos_45 = 45.0f64.to_radians().cos();
                let sin_45 = 45.0f64.to_radians().sin();
                let tan_45 = 45.0f64.to_radians().tan();
                
                let z_rup_bot = s.z_tor + s.w * sin_45;
                
                // Thresholds for R_rup
                // T1: Boundary between Region 1 and 2 (Rx = Z_tor * tan(delta))
                //     At this point, R_rup = Z_tor / cos(delta)
                let t1 = s.z_tor / cos_45;
                
                // T2: Boundary between Region 2 and 3 (Rx = W/cos + Z_tor*tan)
                //     At this point, R_rup = W * tan(delta) + Z_tor / cos(delta)
                let t2 = s.w * tan_45 + s.z_tor / cos_45;

                let r_x = if r_rup <= t1 {
                    (r_rup.powi(2) - s.z_tor.powi(2)).sqrt()
                } else if r_rup <= t2 {
                    (r_rup - s.z_tor * cos_45) / sin_45
                } else {
                    // Region 3
                    // R_rup^2 = (Rx - W_proj)^2 + Z_rup_bot^2
                    // Rx = sqrt(R_rup^2 - Z_rup_bot^2) + W_proj
                    let val = r_rup.powi(2) - z_rup_bot.powi(2);
                    if val < 0.0 {
                        // Should not happen if r_rup > t2
                        0.0 
                    } else {
                        val.sqrt() + s.w * cos_45
                    }
                };
                
                let w_proj = s.w * cos_45;
                let r_jb = if r_x <= w_proj {
                    0.0
                } else {
                    r_x - w_proj
                };

                let eq = EQSource::new(
                    true,
                    s.m,
                    1,
                    0,
                    0.0, 0.0, // Center at 0,0
                    Some(s.w),
                    None, // Length unknown
                    (0.70710678, 0.0, 0.70710678), // Normal for 45 deg dip towards +X.
                    90.0, // Reverse
                    true, // Hanging Wall Effect
                    Some(s.z_hyp),
                    Some(s.z_tor),
                    Some(s.z_bot),
                    Region::California,
                );

                let site = Site::new(
                    0,
                    0.0, 0.0, // Dummy coords
                    0.0,
                    1.0, 
                    760.0,
                    Some(0.6068),
                    false 
                ).with_distances(Some(r_rup), Some(r_jb), Some(r_x));
                
                let pga = cb14.calc_log_y(&eq, &site, 0.0, IMType::PGA).exp();
                row.push_str(&format!(",{:.6}", pga));
            }
            writeln!(file, "{}", row).unwrap();
        }
    }

    #[test]
    fn test_cb14_std() {
        let cb14 = CB14::new();
        
        // Data from Table 2
        // T, sigma (M<=4.5), sigma (M>=5.5)
        let data = vec![
            (0.010, 0.838, 0.590),
            (0.020, 0.848, 0.594),
            (0.030, 0.870, 0.609),
            (0.050, 0.928, 0.642),
            (0.075, 0.930, 0.679),
            (0.10, 0.888, 0.690),
            (0.15, 0.859, 0.667),
            (0.20, 0.833, 0.647),
            (0.25, 0.818, 0.630),
            (0.30, 0.803, 0.642),
            (0.40, 0.776, 0.649),
            (0.50, 0.764, 0.665),
            (0.75, 0.743, 0.712),
            (1.0, 0.746, 0.720),
            (1.5, 0.735, 0.723),
            (2.0, 0.727, 0.711),
            (3.0, 0.726, 0.713),
            (4.0, 0.753, 0.683),
            (5.0, 0.733, 0.693),
            (7.5, 0.695, 0.700),
            (10.0, 0.642, 0.698),
        ];

        // Site condition: Linear (Vs30 >= k1). k1 ranges ~800-1100. 
        // Using Vs30 = 1500.0 to be safe and ensure alpha = 0.
        let site = Site::new(0, 0.0, 0.0, 0.0, 0.0, 1500.0, Some(0.0), false);

        // Case 1: M <= 4.5 (e.g., M=4.0)
        let eq_low_m = EQSource::new(true, 4.0, 0, 0, 0.0, 0.0, Some(0.0), Some(0.0), (0.0, 1.0, 0.0), 0.0, false, Some(0.0), Some(0.0), Some(0.0), Region::Global);

        // Case 2: M >= 5.5 (e.g., M=6.0)
        let eq_high_m = EQSource::new(true, 6.0, 0, 0, 0.0, 0.0, Some(0.0), Some(0.0), (0.0, 1.0, 0.0), 0.0, false, Some(0.0), Some(0.0), Some(0.0), Region::Global);

        for (t, sigma_ref_l, sigma_ref_h) in data {
            let im_type = IMType::PSA(t);
            
            // Low M
            let (sigma_l, _, _) = cb14.calc_std_log_y(&eq_low_m, &site, im_type);
            assert!((sigma_l - sigma_ref_l).abs() < 1e-3, "T={}: sigma_low mismatch. Got {}, expected {}", t, sigma_l, sigma_ref_l);

            // High M
            let (sigma_h, _, _) = cb14.calc_std_log_y(&eq_high_m, &site, im_type);
            assert!((sigma_h - sigma_ref_h).abs() < 1e-3, "T={}: sigma_high mismatch. Got {}, expected {}", t, sigma_h, sigma_ref_h);
        }
        
        // Check PGA and PGV separately
        let pga_data = (0.840, 0.588);
        let pgv_data = (0.728, 0.576);
        
        // PGA
        {
            let (sigma_ref_l, sigma_ref_h) = pga_data;
            let im_type = IMType::PGA;
             // Low M
            let (sigma_l, _, _) = cb14.calc_std_log_y(&eq_low_m, &site, im_type);
            assert!((sigma_l - sigma_ref_l).abs() < 1e-3, "PGA: sigma_low mismatch. Got {}, expected {}", sigma_l, sigma_ref_l);
            
            // High M
            let (sigma_h, _, _) = cb14.calc_std_log_y(&eq_high_m, &site, im_type);
            assert!((sigma_h - sigma_ref_h).abs() < 1e-3, "PGA: sigma_high mismatch. Got {}, expected {}", sigma_h, sigma_ref_h);
        }

        // PGV
        {
            let (sigma_ref_l, sigma_ref_h) = pgv_data;
            let im_type = IMType::PGV;
             // Low M
            let (sigma_l, _, _) = cb14.calc_std_log_y(&eq_low_m, &site, im_type);
            assert!((sigma_l - sigma_ref_l).abs() < 1e-3, "PGV: sigma_low mismatch. Got {}, expected {}", sigma_l, sigma_ref_l);
            
            // High M
            let (sigma_h, _, _) = cb14.calc_std_log_y(&eq_high_m, &site, im_type);
            assert!((sigma_h - sigma_ref_h).abs() < 1e-3, "PGV: sigma_high mismatch. Got {}, expected {}", sigma_h, sigma_ref_h);
        }

    }
}
