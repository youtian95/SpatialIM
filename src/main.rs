
use _spatialim::run_simulation;
use std::env;

/// 程序入口
///
/// # 命令行参数
/// 程序接受三个命令行参数，依次为：
/// 1. **震源信息文件路径** (例如 `EQSource.txt`)
/// 2. **场地信息文件路径** (例如 `SiteFile.txt`)
/// 3. **GMPE模型名称** (可选，默认为 `CB14`)
///
/// # 文件格式说明
///
/// ## 1. 震源信息文件 (`EQSource.txt` 或 `EQSource.json`)
/// 支持两种格式：
///
/// ### 文本格式 (旧版)
/// 逐行包含以下字段（空格分隔）：
/// - `ifmedian`: (0/1) 是否输出中位值
/// - `M`: 震级
/// - `N_sim`: 模拟次数
/// - `seed`: (int) 随机数种子
/// - `lon_0`, `lat_0`: 震中经纬度 (°)
/// - `W`: 断层破裂面的倾斜宽度 (down-dip width)。如果未知输入 999。
/// - `length`: 断裂面水平方向长度
/// - `RuptureNormal_x`, `RuptureNormal_y`, `RuptureNormal_z`: 断裂面朝上的法线方向（向东为x, 向北为y, 向上为z）
/// - `lambda`: 滑动角 (rake angle, degree) - 断裂面上滑动的平均角度
/// - `Fhw`: 上盘效应 (hanging wall effect)
///     - `1`: 包含
///     - `0`: 不包含
/// - `Zhyp`: (km) 震源深度 (海平面以下)。如果未知输入 999。
/// - `region`: 区域代码
///     - `0`: 全球 (含台湾)
///     - `1`: 加利福尼亚
///     - `3`: 中国或土耳其
///     - `4`: 意大利
/// - `nPCs`: IM相关性PCA方法模拟考虑的主成分阶数 (推荐 >= 5)
///
/// ### JSON 格式 (推荐)
/// 示例：
/// ```json
/// {
///   "ifmedian": true,
///   "m": 7.0,
///   "n_sim": 1,
///   "seed": 12345,
///   "lon_0": 100.0,
///   "lat_0": 30.0,
///   "rupture_normal": [0.0, 1.0, 1.0],
///   "lambda": 0.0,
///   "fhw": true,
///   "region": "Global",
///   "n_pcs": 5
///   // 可选参数 (w, length, zhyp, z_tor, z_bot) 如果未知可省略
/// }
/// ```
///
/// ## 2. 场地信息文件 (`SiteFile.txt` 或 `SiteFile.csv`)
/// 支持两种格式：
///
/// ### 文本格式 (空格分隔，无标题)
/// 每行代表一个场地的数据，包含以下字段（空格分隔）：
/// - `ID`: (int) 场地ID
/// - `lon`: 经度
/// - `lat`: 纬度
/// - `elevation_km`: 高程 (km)
/// - `period1`: 基本周期
/// - `Vs30_mpers`: 剪切波速 (m/s)
/// - `Z25_km`: 2.5 km/s 剪切波速层的深度 (km)。如果在加利福尼亚或日本且 Z2.5 未知，输入 999。
///
/// ### CSV 格式 (逗号分隔，带标题)
/// 第一行为标题行，程序会根据标题名称自动识别列（支持乱序）。
/// 支持的标题名称（不区分大小写）：
/// - `ID`
/// - `lon` / `longitude`
/// - `lat` / `latitude`
/// - `elevation_km` / `elevation` / `elev`
/// - `period1` / `period`
/// - `Vs30_mpers` / `Vs30`
/// - `Z25_km` / `Z25`
///
/// 示例：
/// ```csv
/// ID,lon,lat,elevation_km,period1,Vs30_mpers,Z25_km
/// 1,120.1,30.1,0.0,1.0,760.0,999.0
/// ```
fn main() {
    // 获取所有参数的迭代器
    let args: Vec<String> = env::args().collect();

    // 检查参数数量
    if args.len() < 3 {
        eprintln!("错误: 参数数量不正确。");
        eprintln!("用法: {} <EQSource.txt> <SiteFile.txt> [GMPE_Model]", args[0]);
        std::process::exit(1);
    }

    let eq_source_path = &args[1];
    let site_file_path = &args[2];
    let gmpe_model = if args.len() > 3 {
        Some(args[3].as_str())
    } else {
        None
    };

    if let Err(e) = run_simulation(eq_source_path, site_file_path, gmpe_model, None) {
        eprintln!("错误: {}", e);
        std::process::exit(1);
    }
}

