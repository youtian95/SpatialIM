use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::fs;
use std::error::Error;
use super::eq_source::{EQSource, Region};
use super::site::Site;
use nalgebra::DMatrix;

/// 解析震源文件
/// 返回 (EQSource, n_pcs)
pub fn parse_eq_source_file(path: &str) -> Result<(EQSource, usize), Box<dyn Error>> {
    // 尝试解析为 JSON
    if path.ends_with(".json") {
        let json_str = fs::read_to_string(path)?;
        let v: serde_json::Value = serde_json::from_str(&json_str)?;
        
        // 尝试从 JSON 中读取 n_pcs，如果不存在则默认为 6
        let n_pcs = v.get("n_pcs")
            .and_then(|val| val.as_u64())
            .map(|val| val as usize)
            .unwrap_or(6);
            
        let eq = EQSource::from_json(&json_str)?;
        return Ok((eq, n_pcs));
    }
    
    let mut file = File::open(path)?;
    let mut content = String::new();
    file.read_to_string(&mut content)?;
    
    // 如果不是 JSON，则按旧格式解析
    // 读取所有非空行
    let lines: Vec<String> = content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|s| s.to_string())
        .collect();

    // 辅助闭包：解析每一行的第一个字段
    let parse_line = |index: usize| -> f64 {
        lines[index]
            .split_whitespace()
            .next()
            .unwrap_or("0")
            .parse()
            .unwrap_or(0.0)
    };

    // 辅助闭包：解析每一行的第一个字段，如果为 999.0 则返回 None
    let parse_line_opt = |index: usize| -> Option<f64> {
        let val = parse_line(index);
        if val == 999.0 { None } else { Some(val) }
    };

    // 特殊处理经纬度行 (第5行: lon_0 lat_0)
    let coords: Vec<f64> = lines[4]
        .split_whitespace()
        .map(|s| s.parse().unwrap_or(0.0))
        .collect();
    let lon_0 = if coords.len() > 0 { coords[0] } else { 0.0 };
    let lat_0 = if coords.len() > 1 { coords[1] } else { 0.0 };

    // 特殊处理法线方向行 (第8行: RuptureNormal_x _y _z)
    let normals: Vec<f64> = lines[7]
        .split_whitespace()
        .map(|s| s.parse().unwrap_or(0.0))
        .collect();
    let normal_x = if normals.len() > 0 { normals[0] } else { 0.0 };
    let normal_y = if normals.len() > 1 { normals[1] } else { 0.0 };
    let normal_z = if normals.len() > 2 { normals[2] } else { 1.0 };

    // 读取 n_pcs (第13行)，如果存在
    let n_pcs = if lines.len() > 12 {
        lines[12].trim().parse::<usize>().unwrap_or(6)
    } else {
        6
    };

    let eq = EQSource::new(
        parse_line(0) as i32 == 1,
        parse_line(1),
        parse_line(2) as usize,
        parse_line(3) as u64,
        lon_0,
        lat_0,
        parse_line_opt(5),
        parse_line_opt(6),
        (normal_x, normal_y, normal_z),
        parse_line(8),
        parse_line(9) as i32 == 1,
        parse_line_opt(10),
        None,
        None,
        Region::from_i32(parse_line(11) as i32),
    );

    Ok((eq, n_pcs))
}

/// 解析场地文件
/// 根据文件后缀调用不同的解析函数
pub fn parse_site_file(path: &str) -> io::Result<Vec<Site>> {
    if path.to_lowercase().ends_with(".csv") {
        parse_site_file_csv(path)
    } else {
        parse_site_file_txt(path)
    }
}

/// 解析空格分隔的文本场地文件 (无标题)
fn parse_site_file_txt(path: &str) -> io::Result<Vec<Site>> {
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut sites = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 7 {
            sites.push(Site {
                id: parts[0].parse().unwrap_or(0),
                lon: parts[1].parse().unwrap_or(0.0),
                lat: parts[2].parse().unwrap_or(0.0),
                elevation_km: parts[3].parse().unwrap_or(0.0),
                period1: parts[4].parse().unwrap_or(0.0),
                vs30: parts[5].parse().unwrap_or(0.0),
                z25: {
                    let val = parts[6].parse().unwrap_or(999.0);
                    if val == 999.0 { None } else { Some(val) }
                },
                r_rup: None,
                r_jb: None,
                r_x: None,
            });
        }
    }
    Ok(sites)
}

/// 解析 CSV 场地文件 (支持标题映射)
fn parse_site_file_csv(path: &str) -> io::Result<Vec<Site>> {
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();
    let mut sites = Vec::new();

    // 读取第一行
    let first_line = match lines.next() {
        Some(Ok(l)) => l,
        _ => return Ok(sites),
    };

    let trimmed_first = first_line.trim();
    if trimmed_first.is_empty() {
        return Ok(sites);
    }

    let first_parts: Vec<&str> = trimmed_first.split(',').map(|s| s.trim()).collect();
    
    // 判断是否有标题：检查第一个字段是否无法解析为数字
    let has_header = first_parts.first().map_or(false, |s| s.parse::<f64>().is_err());

    // 默认列索引: ID, lon, lat, elevation_km, period1, Vs30_mpers, Z25_km
    let mut col_indices = [0, 1, 2, 3, 4, 5, 6]; 

    if has_header {
        // 建立标题映射
        let header_map: std::collections::HashMap<String, usize> = first_parts
            .iter()
            .enumerate()
            .map(|(i, s)| (s.to_lowercase(), i))
            .collect();

        let find_idx = |names: &[&str], default: usize| -> usize {
            for name in names {
                if let Some(&idx) = header_map.get(*name) {
                    return idx;
                }
            }
            default
        };
        
        col_indices[0] = find_idx(&["id"], 0);
        col_indices[1] = find_idx(&["lon", "longitude"], 1);
        col_indices[2] = find_idx(&["lat", "latitude"], 2);
        col_indices[3] = find_idx(&["elevation_km", "elevation", "elev"], 3);
        col_indices[4] = find_idx(&["period1", "period"], 4);
        col_indices[5] = find_idx(&["vs30_mpers", "vs30"], 5);
        col_indices[6] = find_idx(&["z25_km", "z25"], 6);
    } else {
        // 如果没有标题，第一行也是数据
        if first_parts.len() >= 7 {
             sites.push(parse_site_from_parts(&first_parts, &col_indices));
        }
    }

    // 处理剩余行
    for line_res in lines {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() { continue; }
        
        let parts: Vec<&str> = trimmed.split(',').map(|s| s.trim()).collect();
        if parts.len() >= 7 {
             sites.push(parse_site_from_parts(&parts, &col_indices));
        }
    }

    Ok(sites)
}

/// 解析单行 CSV 字段为 Site
///  - `parts`: 字段切片
///  - `indices`: 各字段在切片中的索引，按顺序为 ID, lon, lat, elevation_km, period1, Vs30_mpers, Z25_km
fn parse_site_from_parts(parts: &[&str], indices: &[usize; 7]) -> Site {
    // 辅助闭包：根据索引获取字段值
    let get_val = |idx: usize| -> f64 {
        if idx < parts.len() {
            parts[idx].parse().unwrap_or(0.0)
        } else {
            0.0
        }
    };
    // 辅助闭包：根据索引获取字段值并解析为整数
    let get_int = |idx: usize| -> i32 {
        if idx < parts.len() {
            parts[idx].parse().unwrap_or(0)
        } else {
            0
        }
    };

    Site {
        id: get_int(indices[0]),
        lon: get_val(indices[1]),
        lat: get_val(indices[2]),
        elevation_km: get_val(indices[3]),
        period1: get_val(indices[4]),
        vs30: get_val(indices[5]),
        z25: {
            let val = get_val(indices[6]);
            if val == 999.0 { None } else { Some(val) }
        },
        r_rup: None,
        r_jb: None,
        r_x: None,
    }
}

/// 保存模拟结果到 CSV 文件
/// 每个周期生成一个文件，包含所有场地在所有模拟次数下的 IM 值
/// 
/// # 参数
/// * `output_dir` - 输出目录
/// * `periods` - 模拟的周期列表
/// * `results` - 模拟结果，一共 n_sim 个元素，每个元素为一个 DMatrix，大小为 [n_sites x n_periods]，行对应场地，列对应周期
/// 
/// # 返回
/// - `Result<(), Box<dyn Error>>` - 成功返回 Ok，失败返回错误信息
/// 
/// # 输出文件格式
/// - 文件名格式: `IM_T{period}.csv`
/// - 内容: 行代表场地，列代表模拟次数
pub fn save_simulation_results(
    output_dir: &str,
    periods: &[f64],
    results: &[DMatrix<f64>],
) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(output_dir)?;
    
    if results.is_empty() {
        return Ok(());
    }

    let n_sims = results.len();
    let n_sites = results[0].nrows();
    let n_periods = periods.len();

    // 验证维度
    if results[0].ncols() != n_periods {
        return Err(format!("结果矩阵列数 ({}) 与周期数 ({}) 不匹配", results[0].ncols(), n_periods).into());
    }

    for (j, &period) in periods.iter().enumerate() {
        // 格式化文件名，例如 IM_T0.01.csv
        let file_path = format!("{}/IM_T{}.csv", output_dir, period);
        let mut file = File::create(file_path)?;

        // 写入表头: Site_ID, Sim_1, Sim_2, ..., Sim_N
        write!(file, "Site_ID")?;
        for k in 0..n_sims {
            write!(file, ",Sim_{}", k + 1)?;
        }
        writeln!(file)?;

        // 写入数据: 每个场地一行
        for i in 0..n_sites {
            write!(file, "{}", i + 1)?; // Site ID (假设从1开始)
            for k in 0..n_sims {
                let val = results[k][(i, j)];
                write!(file, ",{:.6}", val)?;
            }
            writeln!(file)?;
        }
    }

    Ok(())
}

/// 保存各场地在其指定周期下（经插值）得到的模拟结果到单独文件。
///
/// 输出格式：CSV，文件名为 `IM_T_site_specified.csv`
/// 列：`Site_ID,Period,Sim_1,Sim_2,...`
/// - 每一行为一个场地（其周期来自场地文件的 `period1` 字段）
/// - `Sim_k` 为第 k 次模拟在该场地周期处的 IM 值
pub fn save_site_period_results(
    output_dir: &str,
    site_periods: &[f64],
    results: &[Vec<f64>],
) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(output_dir)?;

    // 结果维度：results.len() = n_sims；每个内部 Vec 的长度 = n_sites
    let n_sims = results.len();
    let n_sites = if n_sims > 0 { results[0].len() } else { site_periods.len() };

    // 基本校验
    if n_sites != site_periods.len() {
        return Err("site_periods 长度与结果中的场地数不一致".into());
    }

    let file_path = format!("{}/IM_T_site_specified.csv", output_dir);
    let mut file = File::create(file_path)?;

    // 写表头
    write!(file, "Site_ID,Period")?;
    for k in 0..n_sims {
        write!(file, ",Sim_{}", k + 1)?;
    }
    writeln!(file)?;

    // 写每个场地的数据行
    for i in 0..n_sites {
        write!(file, "{}", i + 1)?; // Site ID（从 1 开始）
        write!(file, ",{:.6}", site_periods[i])?;
        for k in 0..n_sims {
            let val = results[k][i];
            write!(file, ",{:.6}", val)?;
        }
        writeln!(file)?;
    }

    Ok(())
}
