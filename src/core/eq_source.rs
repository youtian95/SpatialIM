use super::geo::{self, Vec3};
use super::site::Site;
use std::f32::consts::PI;
use serde::{Serialize, Deserialize};

/// 区域代码枚举
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Region {
    Global = 0,
    California = 1,
    ChinaTurkey = 3,
    Italy = 4,
    Japan = 5,
}

impl Region {
    pub fn from_i32(v: i32) -> Self {
        match v {
            0 => Region::Global,
            1 => Region::California,
            3 => Region::ChinaTurkey,
            4 => Region::Italy,
            5 => Region::Japan,
            _ => Region::Global,
        }
    }
}

/// 震源参数结构体
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EQSource {
    /// 是否计算中值 (1: 是, 0: 否)
    pub ifmedian: bool,
    /// 矩震级 (Moment Magnitude)
    pub m: f32,
    /// 模拟次数
    pub n_sim: usize,
    /// 随机数种子
    pub seed: u64,
    /// 震中经度 (度)
    pub lon_0: f32,
    /// 震中纬度 (度)
    pub lat_0: f32,
    /// 断层宽度 (km)
    #[serde(default)]
    pub w: Option<f32>,
    /// 断层长度 (km)
    #[serde(default)]
    pub length: Option<f32>,
    /// 断层破裂面法向量 (x, y, z)，Z 方向向上
    #[serde(default)]
    pub rupture_normal: Option<(f32, f32, f32)>,
    /// 走向 (Strike, 度)
    #[serde(default)]
    pub strike: Option<f32>,
    /// 倾角 (Dip, 度)
    #[serde(default)]
    pub dip: Option<f32>,
    /// 滑动角 (Rake Angle, 度)。取值范围 [-180, 180]。
    /// * 0: 左旋走滑 (Left-lateral strike-slip)
    /// * 90: 逆断层 (Reverse)
    /// * -90: 正断层 (Normal)
    /// * 180/-180: 右旋走滑 (Right-lateral strike-slip)
    pub lambda: f32,
    /// 是否考虑上盘效应 (Hanging Wall Effect)
    pub fhw: bool,
    /// 震源深度 (km)
    #[serde(default)]
    pub zhyp: Option<f32>,
    /// 断层顶部深度 (km)，如果未知则为 None
    #[serde(default)]
    pub z_tor: Option<f32>,
    /// 断层底部深度 (km)，如果未知则为 None
    #[serde(default)]
    pub z_bot: Option<f32>,
    /// 区域代码
    pub region: Region,
}

impl EQSource {
    /// 创建新的 EQSource 实例并自动计算几何参数
    ///
    /// # 参数说明
    /// 未知参数请传入`None`，程序将自动估算。 断层宽度 w、 断层长度 length、断层顶部深度 z_tor 可以直接根据震级估算或者直接提供准确值，然后根据宽度 w 计算断层底部深度 z_bot。如果震源深度 zhyp 未知，那么根据 z_tor 和 z_bot 计算震源深度 zhyp。
    ///
    /// * `ifmedian` - 是否计算中值 (必须)
    /// * `m` - 矩震级 (必须)
    /// * `n_sim` - 模拟次数 (必须)
    /// * `seed` - 随机数种子 (必须)
    /// * `lon_0` - 震中经度 (必须)
    /// * `lat_0` - 震中纬度 (必须)
    /// * `w` - 断层宽度 (可选)。如果未知，请传入 `None`，程序将自动估算。
    /// * `length` - 断层长度 (可选)。如果未知，请传入 `None`，程序将自动估算。
    /// * `rupture_normal` - 断层破裂面法向量 (可选)，或者提供 Strike 和 Dip。如果都提供了，以 `rupture_normal` 为准。
    /// * `strike` - 走向 (可选)
    /// * `dip` - 倾角 (可选)
    /// * `lambda` - 滑动角 (必须)
    /// * `fhw` - 是否考虑上盘效应 (必须)
    /// * `zhyp` - 震源深度 (可选)。如果未知，请传入 `None`，程序将自动估算。
    /// * `region` - 区域代码 (必须)
    /// 
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ifmedian: bool,
        m: f32,
        n_sim: usize,
        seed: u64,
        lon_0: f32,
        lat_0: f32,
        w: Option<f32>,
        length: Option<f32>,
        rupture_normal: Option<(f32, f32, f32)>,
        strike: Option<f32>,
        dip: Option<f32>,
        lambda: f32,
        fhw: bool,
        zhyp: Option<f32>,
        z_tor: Option<f32>,
        z_bot: Option<f32>,
        region: Region,
    ) -> Self {
        let mut eq = EQSource {
            ifmedian,
            m,
            n_sim,
            seed,
            lon_0,
            lat_0,
            w,
            length,
            rupture_normal,
            strike,
            dip,
            lambda,
            fhw,
            zhyp,
            z_tor,
            z_bot,
            region,
        };
        
        eq.estimate_unknown_parameters();
        
        eq
    }

    /// 从 JSON 字符串解析 EQSource，并自动估算未知参数
    pub fn from_json(json_str: &str) -> Result<Self, serde_json::Error> {
        let mut eq: Self = serde_json::from_str(json_str)?;
        eq.estimate_unknown_parameters();
        Ok(eq)
    }

    fn estimate_unknown_parameters(&mut self) {
        // 0. Estimate Rupture Normal if missing
        if self.rupture_normal.is_none() {
            if let (Some(strike), Some(dip)) = (self.strike, self.dip) {
                let n = Self::strike_dip_2_rupture_normal(strike, dip);
                self.rupture_normal = Some((n.x, n.y, n.z));
            } else {
                panic!("Must provide either 'rupture_normal' or both 'strike' and 'dip'.");
            }
        }

        let delta = self.delta();
        let sin_delta = delta.to_radians().sin();

        // Internal helper functions

        // Wells and Coppersmith (1994) - Length
        let calc_length = |m: f32, lambda: f32| -> f32 {
            let (a, b) = if lambda.abs() <= 45.0 || lambda.abs() >= 135.0 {
                // strike-slip
                (-3.55, 0.74)
            } else if lambda > 45.0 && lambda < 135.0 {
                // reverse
                (-2.86, 0.63)
            } else if lambda < -45.0 && lambda > -135.0 {
                // normal
                (-2.01, 0.50)
            } else {
                // unknown
                (-3.22, 0.69)
            };
            10.0f32.powf(a + b * m)
        };

        // Wells and Coppersmith (1994) - Width
        let calc_width_wc94 = |m: f32| -> f32 {
            let w_wc = 10.0f32.powf((m - 4.07) / 0.98).sqrt();
             
            w_wc
        };

        // Kaklamanos et al. (2011) - Width
        // J Kaklamanos, L G Baise, D M Boore. Estimating Unknown Input Parameters when Implementing the NGA Ground-Motion Prediction Equations in Engineering Practice. Earthquake Spectra, 2011, 27(4): 1219-1235.

        let _calc_width_kaklamanos = |mag: f32, lambda: f32| -> f32 {
            let abs_lambda = lambda.abs();
            if abs_lambda <= 45.0 || abs_lambda >= 135.0 {
                // strike-slip
                10.0f32.powf(-0.76 + 0.27 * mag)
            } else if lambda > 45.0 && lambda < 135.0 {
                // reverse
                10.0f32.powf(-1.61 + 0.41 * mag)
            } else {
                // normal
                10.0f32.powf(-1.14 + 0.35 * mag)
            }
        };

        // Equations 4 and 5 in Chiou and Youngs (2014)
        let calc_z_tor = |m: f32, is_reverse: bool| -> f32 {
            if is_reverse {
                let val = (m - 5.849).max(0.0);
                (2.704 - 1.226 * val).max(0.0).powi(2)
            } else {
                let val = (m - 4.970).max(0.0);
                (2.673 - 1.136 * val).max(0.0).powi(2)
            }
        };

        // Campbell and Bozorgnia (2013) - Z_HYP
        let calc_zhyp = |m: f32, delta: f32, z_bot: f32, z_tor: f32, _w: f32, _sin_delta: f32| -> f32 {
            let z_bot_eff = z_bot;

            let f_dz_m = if m < 6.75 {
                -4.317 + 0.984 * m
            } else {
                2.325
            };

            let f_dz_delta = if delta <= 40.0 {
                0.0445 * (delta - 40.0)
            } else {
                0.0
            };

            let term1 = f_dz_m + f_dz_delta;
            let term2 = (0.9 * (z_bot_eff - z_tor)).ln();
            
            let ln_dz = term1.min(term2);
            z_tor + ln_dz.exp()
        };

        // Apply estimations

        // 1. Estimate Length
        if self.length.is_none() || self.length == Some(0.0) {
            self.length = Some(calc_length(self.m, self.lambda));
        }

        // 2. Estimate Width
        if self.w.is_none() || self.w == Some(0.0) {
            self.w = Some(calc_width_wc94(self.m));
        }

        // 3. Estimate Z_TOR
        if self.z_tor.is_none() {
            self.z_tor = Some(calc_z_tor(self.m, self.get_f_rv()));
        }

        // 4. Estimate Z_BOT
        if self.z_bot.is_none() {
            // Safe to unwrap because we just estimated them if they were None
            let z_tor = self.z_tor.unwrap();
            let w = self.w.unwrap();
            self.z_bot = Some(z_tor + w * sin_delta);
        }

        // 5. Estimate Z_HYP
        if self.zhyp.is_none() {
            let z_bot = self.z_bot.unwrap();
            let z_tor = self.z_tor.unwrap();
            let w = self.w.unwrap();
            self.zhyp = Some(calc_zhyp(self.m, delta, z_bot, z_tor, w, sin_delta));
        }
    }


    /// 计算断层倾角 (Dip Angle, delta)
    pub fn delta(&self) -> f32 {
        let rn = self.rupture_normal.expect("Rupture normal should be set");
        let normal = Vec3::new(rn.0, rn.1, rn.2).normalize();
        let z_unit = Vec3::new(0.0, 0.0, 1.0);
        (normal.dot(&z_unit).abs()).acos() / PI * 180.0
    }

    /// F_RV: 是否为逆断层 (reverse and reverse-oblique faulting)
    pub fn get_f_rv(&self) -> bool {
        self.lambda > 30.0 && self.lambda < 150.0
    }

    /// F_NM: 是否为正断层 (normal and normal-oblique faulting)
    pub fn get_f_nm(&self) -> bool {
        self.lambda < -30.0 && self.lambda > -150.0
    }

    /// 是否为日本区域
    pub fn get_if_is_japan_site(&self) -> bool {
        self.region == Region::Japan
    }

    /// 计算 Rrup (Rupture Distance): 场地到断层破裂面的最短距离，单位 km
    pub fn calc_rrup(&self, site: &Site) -> f32 {
        let (x, y) = geo::latlon2xy(site.lon, site.lat, self.lon_0, self.lat_0);
        
        // Site Z relative to Hypocenter (which is at 0,0,0 in rupture_4_points)
        // Hypocenter is at depth zhyp. Site is at elevation E.
        // Relative Z = E + zhyp.
        let z = site.elevation_km + self.zhyp.expect("Zhyp should be estimated");
        
        let site_p = Vec3::new(x, y, z);
        
        let rupture_points = self.rupture_4_points();
        let rn = self.rupture_normal.expect("Rupture normal should be set");
        let normal = Vec3::new(rn.0, rn.1, rn.2);
        
        geo::calc_rrup(site_p, &rupture_points, normal)
    }

    /// 计算 Rjb (Joyner-Boore Distance): 场地到断层破裂面在地表投影的最短距离，单位 km
    pub fn calc_rjb(&self, site: &Site) -> f32 {
        let (x, y) = geo::latlon2xy(site.lon, site.lat, self.lon_0, self.lat_0);
        let z = 0.0;
        let site_p = Vec3::new(x, y, z);
        
        let rupture_points = self.rupture_4_points();
        
        geo::calc_rjb(site_p, &rupture_points)
    }

    /// 计算 Rx: 场地到断层迹线（或其延伸线）的水平垂直距离，单位 km
    pub fn calc_rx(&self, site: &Site) -> f32 {
        let (x, y) = geo::latlon2xy(site.lon, site.lat, self.lon_0, self.lat_0);
        let z = 0.0;
        let site_p = Vec3::new(x, y, z);
        
        let rupture_points = self.rupture_4_points();
        let rn = self.rupture_normal.expect("Rupture normal should be set");
        let normal = Vec3::new(rn.0, rn.1, rn.2);
        
        geo::calc_rx(site_p, &rupture_points, normal)
    }

    /// 计算断层破裂面的四个角点坐标，坐标以震源(Hypocenter)为原点，单位 km
    ///
    /// 返回点的顺序构成一个闭合矩形环：
    /// * P1: 上边缘 (Up-Dip) + 走向正向 (Strike Positive)
    /// * P2: 上边缘 (Up-Dip) + 走向负向 (Strike Negative)
    /// * P3: 下边缘 (Down-Dip) + 走向负向 (Strike Negative)
    /// * P4: 下边缘 (Down-Dip) + 走向正向 (Strike Positive)
    ///
    /// 即：P1 -> P2 (上边缘), P3 -> P4 (下边缘)
    pub fn rupture_4_points(&self) -> Vec<Vec3> {
        let rn = self.rupture_normal.expect("Rupture normal should be set");
        let normal = Vec3::new(rn.0, rn.1, rn.2).normalize();
        let z_unit = Vec3::new(0.0, 0.0, 1.0);
        let mut hor_unit = z_unit.cross(&normal);
        
        if hor_unit.norm_squared() < 1e-9 {
            panic!("Horizontal fault (normal parallel to Z axis) is not supported.");
        } else {
            hor_unit = hor_unit.normalize();
        }
        
        // Up-Dip vector (in the plane of the fault, pointing up)
        let ramp_up_dir = normal.cross(&hor_unit).normalize();
        
        // Calculate distances from Hypocenter (0,0,0)
        
        // Along Strike: Assuming Hypocenter is centered
        let half_len = self.length.expect("Length should be estimated") / 2.0;
        
        // Along Dip:
        // Hypocenter depth = zhyp. Top depth = z_tor.
        // Distance from Top Edge to Hypocenter (along dip)
        let sin_delta = self.delta().to_radians().sin();
        
        let dist_top_to_hyp = if sin_delta.abs() > 1e-6 {
            (self.zhyp.expect("Zhyp should be estimated") - self.z_tor.expect("Ztor should be estimated")) / sin_delta
        } else {
            self.w.expect("Width should be estimated") / 2.0 
        };
        
        let dist_hyp_to_bot = self.w.expect("Width should be estimated") - dist_top_to_hyp;
        
        // Vectors from Hypocenter
        let v_strike = hor_unit * half_len;
        let v_up = ramp_up_dir * dist_top_to_hyp;
        let v_down = ramp_up_dir * (-dist_hyp_to_bot);

        // P1 (Top, Strike+)
        let p1 = v_up + v_strike;
        // P2 (Top, Strike-)
        let p2 = v_up - v_strike;
        // P3 (Bot, Strike-)
        let p3 = v_down - v_strike;
        // P4 (Bot, Strike+)
        let p4 = v_down + v_strike;
        
        vec![p1, p2, p3, p4]
    }

    /// 将走向 (Strike) 和 倾角 (Dip) 转换为断层破裂面法向量
    /// 
    /// 输入：
    ///  - Strike: 顺时针从北开始 (0-360)
    ///  - Dip: 从水平面向下 (0-90)
    /// 
    /// 返回
    ///  - 法向量 (x, y, z)，其中 z 向上
    pub fn strike_dip_2_rupture_normal(strike_deg: f32, dip_deg: f32) -> Vec3 {
        let strike_rad = strike_deg.to_radians();
        let dip_rad = dip_deg.to_radians();
        
        // Normal vector (pointing into hanging wall, i.e., Up)
        // n_x = sin(delta) * cos(phi)
        // n_y = -sin(delta) * sin(phi)
        // n_z = cos(delta)
        
        let nx = dip_rad.sin() * strike_rad.cos();
        let ny = -dip_rad.sin() * strike_rad.sin();
        let nz = dip_rad.cos();
        
        Vec3::new(nx, ny, nz)
    }

    /// 将断层破裂面法向量转换为走向 (Strike) 和 倾角 (Dip)
    /// 返回: (Strike, Dip) 单位为度
    pub fn rupture_normal_2_strike_dip(normal: Vec3) -> (f32, f32) {
        let n = normal.normalize();
        // 确保法向量指向上方 (z >= 0)
        let n = if n.z < -1e-6 { -n } else { n };
        
        let dip_rad = n.z.acos();
        let dip_deg = dip_rad.to_degrees();
        
        // 如果倾角接近 0，走向未定义 (返回 0)
        if dip_deg.abs() < 1e-6 {
            return (0.0, 0.0);
        }
        
        // nx = sin(delta) * cos(phi)
        // ny = -sin(delta) * sin(phi)
        // cos(phi) = nx / sin(delta)
        // sin(phi) = -ny / sin(delta)
        // phi = atan2(-ny, nx)
        
        let sin_dip = dip_rad.sin();
        let phi_rad = (-n.y / sin_dip).atan2(n.x / sin_dip);
        let mut strike_deg = phi_rad.to_degrees();
        
        if strike_deg < 0.0 {
            strike_deg += 360.0;
        }
        
        (strike_deg, dip_deg)
    }


}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn create_dummy_eq_source() -> EQSource {
        EQSource::new(
            true,
            7.0,
            1,
            12345,
            100.0,
            30.0,
            Some(10.0),
            Some(20.0),
            Some((0.0, 1.0, 1.0)), // 45 degree dip
            None,
            None,
            0.0,
            true,
            Some(10.0),
            None,
            None,
            Region::Global,
        )
    }

    fn create_dummy_site(lon: f32, lat: f32) -> Site {
        Site {
            id: 1,
            lon,
            lat,
            elevation_km: 0.0,
            period1: 0.0,
            vs30: 760.0,
            z25: Some(1.0),
            r_rup: None,
            r_jb: None,
            r_x: None,
        }
    }

    #[test]
    fn test_delta() {
        let mut eq = create_dummy_eq_source();
        
        // Vertical fault (normal along Y)
        eq.rupture_normal = Some((0.0, 1.0, 0.0));
        assert_relative_eq!(eq.delta(), 90.0, epsilon = 1e-6);

        // Horizontal fault (normal along Z)
        eq.rupture_normal = Some((0.0, 0.0, 1.0));
        assert_relative_eq!(eq.delta(), 0.0, epsilon = 1e-6);

        // 45 degree dip
        eq.rupture_normal = Some((0.0, 1.0, 1.0));
        assert_relative_eq!(eq.delta(), 45.0, epsilon = 1e-6);
    }

    #[test]
    fn test_rupture_4_points() {
        let mut eq = create_dummy_eq_source();
        eq.length = Some(20.0);
        eq.w = Some(10.0);
        // fault along X axis
        eq.rupture_normal = Some((0.0, 1.0, 1.0));
        
        let points = eq.rupture_4_points();
        assert_eq!(points.len(), 4);

        // Check dimensions
        let p1 = points[0];
        let p2 = points[1];
        let p3 = points[2];
        let p4 = points[3];

        // P1 -> P2 is length (along X)
        assert_relative_eq!((p1 - p2).norm(), 20.0, epsilon = 1e-6);
        // P2 -> P3 is width (along Z, since it's vertical)
        assert_relative_eq!((p2 - p3).norm(), 10.0, epsilon = 1e-6);
        // P3 -> P4 is length
        assert_relative_eq!((p3 - p4).norm(), 20.0, epsilon = 1e-6);
        // P4 -> P1 is width
        assert_relative_eq!((p4 - p1).norm(), 10.0, epsilon = 1e-6);
    }

    #[test]
    fn test_calc_rrup() {
        let mut eq = create_dummy_eq_source();
        eq.zhyp = Some(10.0);
        eq.rupture_normal = Some((0.0, 1.0, 0.0)); // Vertical fault
        eq.w = Some(10.0); 
        eq.length = Some(10.0); 
        // Set z_tor to make fault centered at zhyp
        eq.z_tor = Some(eq.zhyp.unwrap() - eq.w.unwrap() / 2.0); // 10 - 5 = 5.0
        
        let site = create_dummy_site(100.0, 30.0); // Site at epicenter
        
        // Site is at (0, 0, 10). (Elevation 0 + zhyp 10)
        // Fault center (Hypocenter) at (0, 0, 0) in local coords.
        // Vertical fault (normal Y) extends in X and Z.
        // Width 10 means Z from -5 to 5 relative to Hypocenter.
        // Site Z=10 is 5km away from top edge Z=5.
        let rrup = eq.calc_rrup(&site);
        assert_relative_eq!(rrup, 5.0, epsilon = 1e-6);
    }

    #[test]
    fn test_calc_rjb() {
        let mut eq = create_dummy_eq_source();
        eq.zhyp = Some(10.0);
        eq.rupture_normal = Some((0.0, 1.0, 0.0)); // Vertical fault
        eq.w = Some(10.0);
        eq.length = Some(10.0);
        eq.z_tor = Some(eq.zhyp.unwrap() - eq.w.unwrap() / 2.0);
        
        // Site at epicenter
        let site = create_dummy_site(100.0, 30.0);
        
        // Rjb should be 0 because site (0,0) is on the projection of the fault (line along X)
        let rjb = eq.calc_rjb(&site);
        assert_relative_eq!(rjb, 0.0, epsilon = 1e-6);
        
        // Site far away
        // 1 degree lat is approx 111km
        let site_far = create_dummy_site(100.0, 31.0); 
        let rjb_far = eq.calc_rjb(&site_far);
        assert!(rjb_far > 100.0);
    }

    #[test]
    fn test_calc_rx() {
        let mut eq = create_dummy_eq_source();
        eq.zhyp = Some(10.0);
        eq.rupture_normal = Some((0.0, 1.0, 1.0)); // 45 degree dip
        eq.w = Some(10.0);
        // Set z_tor consistent with centered hypocenter
        let sin_delta = 45.0f32.to_radians().sin();
        eq.z_tor = Some(eq.zhyp.unwrap() - (eq.w.unwrap() / 2.0) * sin_delta);
        
        // Check dip direction
        // Normal (0, 1, 1). Z (0, 0, 1).
        // Strike = Z x N = (-1, 0, 0).
        // Up-Dip = N x Strike = (0, 1, 1) x (-1, 0, 0) = (0, -1, 1).
        // So Up-Dip is towards -Y.
        // Down-Dip is towards +Y.
        // Hanging wall is on +Y side.
        
        let site_hw = create_dummy_site(100.0, 30.1); // Site at +Y (Hanging Wall)
        let rx_hw = eq.calc_rx(&site_hw);
        
        // Calculate expected Rx
        // Site Y relative to epicenter
        let (_, site_y) = geo::latlon2xy(100.0, 30.1, 100.0, 30.0);
        // Fault top edge Y relative to epicenter
        // Width = 10. Half width = 5.
        // Dip 45 deg. Horizontal projection of half width = 5 * cos(45).
        // Normal is (0, 1, 1), so dip direction is +Y.
        // Up-dip direction is -Y.
        // Top edge is at -5 * cos(45).
        let fault_trace_y = -5.0 * (45.0_f32.to_radians()).cos();
        
        let expected_rx = site_y - fault_trace_y;
        assert_relative_eq!(rx_hw, expected_rx, epsilon = 1e-5);
        
        let site_fw = create_dummy_site(100.0, 29.9); // Site at -Y (Foot Wall)
        let rx_fw = eq.calc_rx(&site_fw);
        
        let (_, site_y_fw) = geo::latlon2xy(100.0, 29.9, 100.0, 30.0);
        let expected_rx_fw = site_y_fw - fault_trace_y;
        assert_relative_eq!(rx_fw, expected_rx_fw, epsilon = 1e-5);
    }

    #[test]
    fn test_geometry_conversions() {
        let epsilon = 1e-6;

        // Test 1: Strike 0 (North), Dip 90 (Vertical)
        // Normal should be East (1, 0, 0)
        let n = EQSource::strike_dip_2_rupture_normal(0.0, 90.0);
        assert_relative_eq!(n.x, 1.0, epsilon = epsilon);
        assert_relative_eq!(n.y, 0.0, epsilon = epsilon);
        assert_relative_eq!(n.z, 0.0, epsilon = epsilon);
        
        let (s, d) = EQSource::rupture_normal_2_strike_dip(n);
        assert_relative_eq!(s, 0.0, epsilon = epsilon);
        assert_relative_eq!(d, 90.0, epsilon = epsilon);

        // Test 2: Strike 90 (East), Dip 90 (Vertical)
        // Normal should be South (0, -1, 0)
        // n_x = sin(90)*cos(90) = 0
        // n_y = -sin(90)*sin(90) = -1
        // n_z = cos(90) = 0
        let n = EQSource::strike_dip_2_rupture_normal(90.0, 90.0);
        assert_relative_eq!(n.x, 0.0, epsilon = epsilon);
        assert_relative_eq!(n.y, -1.0, epsilon = epsilon);
        assert_relative_eq!(n.z, 0.0, epsilon = epsilon);

        let (s, d) = EQSource::rupture_normal_2_strike_dip(n);
        assert_relative_eq!(s, 90.0, epsilon = epsilon);
        assert_relative_eq!(d, 90.0, epsilon = epsilon);

        // Test 3: Strike 0, Dip 0 (Horizontal)
        // Normal should be Up (0, 0, 1)
        let n = EQSource::strike_dip_2_rupture_normal(0.0, 0.0);
        assert_relative_eq!(n.x, 0.0, epsilon = epsilon);
        assert_relative_eq!(n.y, 0.0, epsilon = epsilon);
        assert_relative_eq!(n.z, 1.0, epsilon = epsilon);

        let (_s, d) = EQSource::rupture_normal_2_strike_dip(n);
        // Strike is undefined for horizontal plane, but dip should be 0
        assert_relative_eq!(d, 0.0, epsilon = epsilon);
    }
}
