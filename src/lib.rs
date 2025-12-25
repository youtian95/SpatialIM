//! 库的主模块


use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError; 


pub mod core;
use core::io::{parse_eq_source_file, parse_site_file};
use core::simulator::Simulator;
use std::error::Error;

/// # 运行模拟的主函数
/// 
/// ## 参数
/// * `eq_source_path` - 震源文件路径
/// * `site_file_path` - 场地文件路径
/// * `gmpe_model` - GMPE 模型名称 (可选，默认为 "CB14")
pub fn run_simulation(eq_source_path: &str, site_file_path: &str, gmpe_model: Option<&str>) -> Result<(), Box<dyn Error>> {
    let model = gmpe_model.unwrap_or("CB14");
    println!("正在处理...");
    println!("  震源文件: {}", eq_source_path);
    println!("  场地文件: {}", site_file_path);
    println!("  GMPE模型: {}", model);

    // 读取并解析文件
    let (eq_source, n_pcs) = parse_eq_source_file(eq_source_path)?;
    let sites = parse_site_file(site_file_path)?;

    // 初始化并运行模拟器
    let sim = Simulator::new(eq_source, sites, model.to_string(), n_pcs);
    sim.run();

    Ok(())
}



// ----------------------------------------------------------------
// Python 接口层 (包装函数)
// ----------------------------------------------------------------

/// 这是专门给 Python 调用的包装函数
/// 它调用上面的纯 Rust 函数，并将 Rust 错误转换为 Python 异常
#[pyfunction]
#[pyo3(name = "run_simulation")] // 在 Python 中显示的名字仍为 run_simulation
fn run_simulation_py(eq_source_path: &str, site_file_path: &str, gmpe_model: Option<&str>) -> PyResult<()> {
    // 调用纯 Rust 函数，并映射错误
    run_simulation(eq_source_path, site_file_path, gmpe_model)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

// ----------------------------------------------------------------
// Python 模块定义
// ----------------------------------------------------------------

#[pymodule]
fn _spatialim(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // 注意这里注册的是包装函数 run_simulation_py
    m.add_wrapped(wrap_pyfunction!(run_simulation_py))?;
    Ok(())
}