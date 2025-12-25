use super::{GMPEModel, GMPEResult};
use crate::core::eq_source::EQSource;
use crate::core::site::Site;

/// 示例：另一个模型 ASK14
pub struct ASK14;

impl ASK14 {
    pub fn new() -> Self {
        ASK14
    }
}

impl GMPEModel for ASK14 {
    fn calc(&self, _eq: &EQSource, _site: &Site) -> Result<GMPEResult, String> {
        println!("Calculating using ASK14 GMPE model (Placeholder)...");
        Ok(GMPEResult {
            psa_median: 0.1,
            psa_sigma: 0.1,
            psa_tau: 0.1,
            psa_phi: 0.1,
            pga_median: 0.1,
            pga_sigma: 0.1,
            pga_tau: 0.1,
            pga_phi: 0.1,
            pgv_median: 0.1,
            pgv_sigma: 0.1,
            pgv_tau: 0.1,
            pgv_phi: 0.1,
        })
    }
}
