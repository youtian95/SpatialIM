//! 库的主模块


use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError; 


pub mod core;
use core::io::{parse_eq_source_file, parse_site_file};
use core::simulator::Simulator;
use core::grid::generate_grid_points_from_bounds;
use core::site::Site;
use core::geo;
use std::error::Error;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

/// # 运行模拟的主函数
/// 
/// ## 参数
/// * `eq_source_path` - 震源文件路径
/// * `site_file_path` - 场地文件路径
/// * `gmpe_model` - GMPE 模型名称 (可选，默认为 "CB14")
/// * `output_dir` - 输出目录 (可选，默认为 "output")
pub fn run_simulation(eq_source_path: &str, site_file_path: &str, gmpe_model: Option<&str>, output_dir: Option<&str>) -> Result<(), Box<dyn Error>> {
    let model = gmpe_model.unwrap_or("CB14");
    let out_dir = output_dir.unwrap_or("output");
    println!("正在处理(直接输入场地模式，如果场地数量太多，会自动生成网格并插值计算)...");
    println!("  震源文件: {}", eq_source_path);
    println!("  场地文件: {}", site_file_path);
    println!("  GMPE模型: {}", model);
    println!("  输出目录: {}", out_dir);

    // 读取并解析文件
    let (eq_source, n_pcs) = parse_eq_source_file(eq_source_path)?;
    let sites = parse_site_file(site_file_path)?;

    // 初始化并运行模拟器
    let sim = Simulator::new(eq_source, sites, model.to_string(), n_pcs)
        .with_output_dir(out_dir.to_string());
    sim.run();

    Ok(())
}

/// # 运行模拟的主函数 (网格模式)
/// 
/// ## 参数
/// * `eq_source_path` - 震源文件路径
/// * `min_lon`, `max_lon`, `min_lat`, `max_lat` - 区域范围
/// * `grid_spacing_km` - 网格间距 (km)
/// * `gmpe_model` - GMPE 模型名称 (可选，默认为 "CB14")
/// * `output_dir` - 输出目录 (可选，默认为 "output")
#[allow(clippy::too_many_arguments)]
pub fn run_simulation_grid(
    eq_source_path: &str,
    min_lon: f32, max_lon: f32, min_lat: f32, max_lat: f32,
    grid_spacing_km: Option<f32>,
    gmpe_model: Option<&str>,
    output_dir: Option<&str>
) -> Result<(), Box<dyn Error>> {
    let model = gmpe_model.unwrap_or("CB14");
    let out_dir = output_dir.unwrap_or("output");
    println!("正在处理(根据区域范围自动划分网格)...");
    println!("  震源文件: {}", eq_source_path);
    println!("  范围: [{}, {}] x [{}, {}]", min_lon, max_lon, min_lat, max_lat);
    println!("  网格间距: {:?} km", grid_spacing_km);
    println!("  GMPE模型: {}", model);
    println!("  输出目录: {}", out_dir);

    // 读取并解析文件
    let (eq_source, n_pcs) = parse_eq_source_file(eq_source_path)?;
    
    // 生成网格
    let (points, _nx, _ny) = generate_grid_points_from_bounds(min_lon, max_lon, min_lat, max_lat, grid_spacing_km);
    println!("  生成网格点数: {}", points.len());

    // 确保输出目录存在
    if !Path::new(out_dir).exists() {
        fs::create_dir_all(out_dir)?;
    }

    // 保存网格坐标信息
    let coords_file_path = Path::new(out_dir).join("grid_coordinates.csv");
    let mut coords_file = File::create(&coords_file_path)?;
    // 写入表头: ID, Longitude, Latitude, LocalX, LocalY
    // LocalX, LocalY 是相对于震源参考点 (eq_source.lon_0, eq_source.lat_0) 的坐标
    writeln!(coords_file, "ID,Longitude,Latitude,LocalX,LocalY")?;

    let lon_0 = eq_source.lon_0;
    let lat_0 = eq_source.lat_0;

    // 创建 Site 对象并写入坐标文件
    let sites: Vec<Site> = points.iter().enumerate().map(|(i, (lon, lat))| {
        Site::new(
            i as i32,
            *lon,
            *lat,
            0.0, // elevation
            0.0, // period1
            760.0, // vs30
            None, // z25
            false // is_japan
        )
    }).collect();

    // 遍历 sites 写入坐标文件
    for site in &sites {
        let (x, y) = geo::latlon2xy(site.lon, site.lat, lon_0, lat_0);
        writeln!(coords_file, "{},{:.6},{:.6},{:.6},{:.6}", site.id, site.lon, site.lat, x, y)?;
    }
    println!("  网格坐标已保存至: {:?}", coords_file_path);

    // 初始化并运行模拟器
    let sim = Simulator::new(eq_source, sites, model.to_string(), n_pcs)
        .with_grid_threshold(usize::MAX) // 网格模式不进行阈值裁剪
        .with_output_dir(out_dir.to_string());
    sim.run();

    Ok(())
}



// ----------------------------------------------------------------
// Python 接口层 (包装函数)
// ----------------------------------------------------------------

/// 这是专门给 Python 调用的包装函数
/// 它调用上面的纯 Rust 函数，并将 Rust 错误转换为 Python 异常
#[pyfunction]
#[pyo3(name = "run_simulation", signature = (eq_source_path, site_file_path, gmpe_model=None, output_dir=None))]
fn run_simulation_py(eq_source_path: &str, site_file_path: &str, gmpe_model: Option<&str>, output_dir: Option<&str>) -> PyResult<()> {
    // 调用纯 Rust 函数，并映射错误
    run_simulation(eq_source_path, site_file_path, gmpe_model, output_dir)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

#[pyfunction]
#[pyo3(name = "run_simulation_grid", signature = (eq_source_path, min_lon, max_lon, min_lat, max_lat, grid_spacing_km=None, gmpe_model=None, output_dir=None))]
fn run_simulation_grid_py(
    eq_source_path: &str,
    min_lon: f32, max_lon: f32, min_lat: f32, max_lat: f32,
    grid_spacing_km: Option<f32>,
    gmpe_model: Option<&str>,
    output_dir: Option<&str>
) -> PyResult<()> {
    run_simulation_grid(eq_source_path, min_lon, max_lon, min_lat, max_lat, grid_spacing_km, gmpe_model, output_dir)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}

// ----------------------------------------------------------------
// Python 模块定义
// ----------------------------------------------------------------

#[pymodule]
fn _spatialim(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // 注意这里注册的是包装函数 run_simulation_py
    m.add_wrapped(wrap_pyfunction!(run_simulation_py))?;
    m.add_wrapped(wrap_pyfunction!(run_simulation_grid_py))?;
    Ok(())
}