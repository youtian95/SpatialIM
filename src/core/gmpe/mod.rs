pub mod cb14;
pub mod ask14;

use self::cb14::CB14;
use self::ask14::ASK14;
use crate::core::eq_source::EQSource;
use crate::core::site::Site;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IMType {
    PGA,
    PGV,
    /// f32 表示周期 T
    PSA(f32),
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct GMPEResult {
    pub psa_median: f32,
    /// PSA 标准差
    pub psa_sigma: f32,
    /// PSA 事件间标准差
    pub psa_tau: f32,
    /// PSA 事件内标准差
    pub psa_phi: f32,
    pub pga_median: f32,
    pub pga_sigma: f32,
    pub pga_tau: f32,
    pub pga_phi: f32,
    pub pgv_median: f32,
    pub pgv_sigma: f32,
    pub pgv_tau: f32,
    pub pgv_phi: f32,
}

/// GMPE 模型特征（接口）
/// 所有具体的 GMPE 模型（如 CB14, ASK14 等）都必须实现这个 trait
pub trait GMPEModel {
    /// 计算给定震源和场地的地震动强度中值和标准差
    fn calc(&self, eq: &EQSource, site: &Site) -> Result<GMPEResult, String>;
}

/// 工厂函数：根据名称创建对应的 GMPE 模型
pub fn create_gmpe_model(name: &str) -> Option<Box<dyn GMPEModel>> {
    match name {
        "CB14" => Some(Box::new(CB14::new())),
        "ASK14" => Some(Box::new(ASK14::new())),
        // 将来在这里添加更多模型...
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::eq_source::Region;

    pub fn create_dummy_eq_source() -> EQSource {
        EQSource::new(
            true,
            7.0,
            1,
            12345,
            0.0,
            0.0,
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

    pub fn create_dummy_site() -> Site {
        Site {
            id: 1,
            lon: 0.1,
            lat: 0.1,
            elevation_km: 0.0,
            period1: 1.0, // T = 1.0s
            vs30: 760.0,
            z25: Some(1.0),
            r_rup: None,
            r_jb: None,
            r_x: None,
        }
    }

    fn verify_model_basic(model: &dyn GMPEModel) {
        let eq = create_dummy_eq_source();
        let site = create_dummy_site();

        let result = model.calc(&eq, &site);
        assert!(result.is_ok(), "Model calculation failed");
        let res = result.unwrap();
        
        println!("PSA Median: {}, Sigma: {}, Tau: {}", res.psa_median, res.psa_sigma, res.psa_tau);
        
        assert!(res.psa_median > 0.0, "Median should be positive");
        assert!(res.psa_sigma > 0.0, "Sigma should be positive");
        assert!(res.psa_tau > 0.0, "Tau should be positive");
    }

    #[test]
    fn test_cb14_basic() {
        let model = create_gmpe_model("CB14").expect("CB14 model not found");
        verify_model_basic(model.as_ref());
    }

    #[test]
    fn test_ask14_basic() {
        let model = create_gmpe_model("ASK14").expect("ASK14 model not found");
        verify_model_basic(model.as_ref());
    }
}
