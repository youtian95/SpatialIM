/// 场地参数结构体
#[derive(Debug, Clone)]
pub struct Site {
    pub id: i32,
    pub lon: f32,
    pub lat: f32,
    pub elevation_km: f32,
    pub period1: f32,
    /// 剪切波速，单位为m/s
    pub vs30: f32,
    /// 深度到剪切波速为2.5km/s的深度，单位为km
    pub z25: Option<f32>,
    /// Rupture distance (km). If None, calculated from geometry.
    pub r_rup: Option<f32>,
    /// Joyner-Boore distance (km). If None, calculated from geometry.
    pub r_jb: Option<f32>,
    /// Horizontal distance from trace (km). If None, calculated from geometry.
    pub r_x: Option<f32>,
}

impl Site {
    /// 创建新的 Site 实例
    ///
    /// 如果 z25 为 None (未知)，将根据 vs30 和 is_japan 自动估算。
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: i32,
        lon: f32,
        lat: f32,
        elevation_km: f32,
        period1: f32,
        vs30: f32,
        z25: Option<f32>,
        is_japan: bool,
    ) -> Self {
        let mut site = Site {
            id,
            lon,
            lat,
            elevation_km,
            period1,
            vs30,
            z25,
            r_rup: None,
            r_jb: None,
            r_x: None,
        };
        
        if site.z25.is_none() {
            site.z25 = Some(site.estimate_z25(is_japan));
        }
        site
    }

    /// Set manual distances
    pub fn with_distances(mut self, r_rup: Option<f32>, r_jb: Option<f32>, r_x: Option<f32>) -> Self {
        self.r_rup = r_rup;
        self.r_jb = r_jb;
        self.r_x = r_x;
        self
    }

    /// 获取 Z2.5 (km)
    ///
    /// 如果 z25 为 None，则根据 Campbell and Bozorgnia (2014) 的经验公式，利用 Vs30 估算 Z2.5。
    pub fn get_z25(&self, is_japan: bool) -> f32 {
        self.z25.unwrap_or_else(|| self.estimate_z25(is_japan))
    }

    /// 估算 Z2.5 (km)
    ///
    /// 根据 Campbell and Bozorgnia (2014) 的经验公式，利用 Vs30 估算 Z2.5。
    ///
    /// # 参数
    /// * `is_japan` - 是否为日本区域。如果是，使用日本的经验公式；否则使用加州的经验公式（通常也适用于其他活跃构造区）。
    pub fn estimate_z25(&self, is_japan: bool) -> f32 {
        let ln_vs30 = self.vs30.ln();
        
        let ln_z25 = if is_japan {
            // Japan: ln Z2.5 = 5.359 - 1.102 ln Vs30
            5.359 - 1.102 * ln_vs30
        } else {
            // California: ln Z2.5 = 7.089 - 1.144 ln Vs30
            7.089 - 1.144 * ln_vs30
        };
        
        ln_z25.exp()
    }
}
