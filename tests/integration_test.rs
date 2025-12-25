use _spatial_im::core::io::{parse_eq_source_file, parse_site_file};
use _spatial_im::core::simulator::Simulator;
use _spatial_im::run_simulation;
use std::path::PathBuf;

#[test]
fn test_main_function_logic() {
    // 获取测试数据文件的路径
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests/fixtures");

    let eq_path = d.join("EQSource.txt");
    let site_path = d.join("SiteFile.txt");

    // 直接调用封装后的 main 逻辑
    let result = run_simulation(
        eq_path.to_str().unwrap(),
        site_path.to_str().unwrap(),
        Some("CB14")
    );

    assert!(result.is_ok(), "Simulation failed: {:?}", result.err());
}

#[test]
fn test_txt_input_flow() {
    // 获取测试数据文件的路径
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests/fixtures");

    let eq_path = d.join("EQSource.txt");
    let site_path = d.join("SiteFile.txt");

    // 1. 解析输入文件
    let (eq_source, n_pcs) = parse_eq_source_file(eq_path.to_str().unwrap())
        .expect("Failed to parse EQSource.txt");
    
    let sites = parse_site_file(site_path.to_str().unwrap())
        .expect("Failed to parse SiteFile.txt");

    // 验证解析结果是否正确
    assert_eq!(eq_source.m, 7.0);
    assert_eq!(eq_source.n_sim, 5);
    assert_eq!(sites.len(), 2);
    assert_eq!(sites[0].id, 1);

    // 2. 初始化模拟器
    let sim = Simulator::new(eq_source, sites, "CB14".to_string(), n_pcs);

    // 3. 运行模拟 (确保不 panic)
    sim.run();
}

#[test]
fn test_json_input_flow() {
    // 获取测试数据文件的路径
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests/fixtures");

    let eq_path = d.join("EQSource.json");
    let site_path = d.join("SiteFile.txt");

    // 1. 解析输入文件 (JSON 格式)
    let (eq_source, n_pcs) = parse_eq_source_file(eq_path.to_str().unwrap())
        .expect("Failed to parse EQSource.json");
    
    let sites = parse_site_file(site_path.to_str().unwrap())
        .expect("Failed to parse SiteFile.txt");

    // 验证解析结果是否正确
    assert_eq!(eq_source.m, 7.0);
    assert_eq!(eq_source.n_sim, 5);
    // 验证自动估算的参数不为 None (因为我们在 JSON 中省略了它们)
    assert!(eq_source.w.is_some());
    assert!(eq_source.length.is_some());
    assert!(eq_source.zhyp.is_some());
    assert!(eq_source.z_tor.is_some());
    assert!(eq_source.z_bot.is_some());

    // 2. 初始化模拟器
    let sim = Simulator::new(eq_source, sites, "CB14".to_string(), n_pcs);

    // 3. 运行模拟
    sim.run();
}

#[test]
fn test_csv_site_file() {
    // 获取测试数据文件的路径
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests/fixtures");

    let site_path = d.join("SiteFile.csv");

    // 解析 CSV 格式的场地文件
    let sites = parse_site_file(site_path.to_str().unwrap())
        .expect("Failed to parse SiteFile.csv");

    // 验证解析结果
    assert_eq!(sites.len(), 2);
    assert_eq!(sites[0].id, 1);
    assert_eq!(sites[0].vs30, 760.0);
    assert_eq!(sites[1].id, 2);
    assert_eq!(sites[1].vs30, 400.0);
}

#[test]
fn test_csv_site_file_shuffled() {
    // 获取测试数据文件的路径
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests/fixtures");

    let site_path = d.join("SiteFile_shuffled.csv");

    // 解析 CSV 格式的场地文件 (乱序标题)
    let sites = parse_site_file(site_path.to_str().unwrap())
        .expect("Failed to parse SiteFile_shuffled.csv");

    // 验证解析结果
    assert_eq!(sites.len(), 2);
    
    // 验证 ID=1 的数据 (lat=30.1, lon=120.1)
    let s1 = sites.iter().find(|s| s.id == 1).unwrap();
    assert_eq!(s1.lat, 30.1);
    assert_eq!(s1.lon, 120.1);
    assert_eq!(s1.vs30, 760.0);
    assert!(s1.z25.is_none());

    // 验证 ID=2 的数据 (lat=30.2, lon=120.2)
    let s2 = sites.iter().find(|s| s.id == 2).unwrap();
    assert_eq!(s2.lat, 30.2);
    assert_eq!(s2.lon, 120.2);
    assert_eq!(s2.vs30, 400.0);
    assert_eq!(s2.z25, Some(1.5));
}

#[test]
fn test_grid_simulation_logic() {
    use _spatial_im::core::site::Site;
    use _spatial_im::core::eq_source::{EQSource, Region};
    use _spatial_im::core::simulator::Simulator;
    use std::path::Path;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    // Helper to read a specific value from CSV
    fn read_csv_val(file_path: &Path, site_idx: usize, sim_idx: usize) -> f64 {
        let file = File::open(file_path).expect("Failed to open CSV");
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        
        // Skip header
        lines.next();
        
        // Skip to site row
        for _ in 0..site_idx {
            lines.next();
        }
        
        let line = match lines.next() {
            Some(Ok(l)) => l,
            _ => panic!("Site row {} not found in file {:?}", site_idx, file_path),
        };
        let parts: Vec<&str> = line.split(',').collect();
        // parts[0] is Site_ID, parts[1] is Sim_1...
        let val_str = parts[sim_idx + 1];
        val_str.parse::<f64>().expect("Parse error")
    }

    // 1. Create dummy EQSource with ifmedian=true (Direct Median)
    let eq_source = EQSource::new(
        true, // ifmedian = true -> Output Median only
        7.0, // m
        1, // n_sim
        12345, // seed
        120.0, // lon_0
        30.0, // lat_0
        Some(10.0), // w
        Some(20.0), // length
        (1.0, 0.0, 0.0), // rupture_normal (dummy)
        0.0, // lambda (rake)
        false, // fhw
        Some(10.0), // zhyp
        Some(0.0), // z_tor
        Some(15.0), // z_bot
        Region::Global, // region
    );

    // 2. Create dummy Sites
    let sites = vec![
        Site::new(1, 120.0, 30.0, 0.0, 0.0, 760.0, None, false),
        Site::new(2, 120.1, 30.0, 0.0, 0.0, 760.0, None, false),
        Site::new(3, 120.0, 30.1, 0.0, 0.0, 760.0, None, false),
        Site::new(4, 120.1, 30.1, 0.0, 0.0, 760.0, None, false),
    ];

    let output_dir = Path::new("output");
    let target_file = output_dir.join("IM_T1.csv");

    // --- Run 1: Direct Simulation ---
    println!("Running Direct Simulation...");
    let sim_direct = Simulator::new(eq_source.clone(), sites.clone(), "CB14".to_string(), 1)
        .with_grid_threshold(1000); // Force direct mode
    sim_direct.run();
    
    // Read all direct values
    let mut direct_vals = Vec::new();
    for i in 0..sites.len() {
        direct_vals.push(read_csv_val(&target_file, i, 0));
    }

    // --- Run 2: Grid Simulation ---
    println!("Running Grid Simulation...");
    let sim_grid = Simulator::new(eq_source.clone(), sites.clone(), "CB14".to_string(), 1)
        .with_grid_threshold(0) // Force grid mode
        .with_grid_spacing(0.5); 
    sim_grid.run();

    // --- Compare ---
    for i in 0..sites.len() {
        let val_grid = read_csv_val(&target_file, i, 0);
        let val_direct = direct_vals[i];
        
        println!("Site {}: Direct={}, Grid={}", i+1, val_direct, val_grid);

        let diff = (val_direct - val_grid).abs();
        let avg = (val_direct + val_grid) / 2.0;
        let rel_err = if avg > 1e-9 { diff / avg } else { 0.0 };
        
        println!("  Difference: {}, Rel Error: {}", diff, rel_err);

        assert!(rel_err < 0.02, "Grid interpolation result for Site {} should be close to direct calculation (within 2%)", i+1);
    }
}

#[test]
fn test_random_median_vs_direct_median() {
    use _spatial_im::core::site::Site;
    use _spatial_im::core::eq_source::{EQSource, Region};
    use _spatial_im::core::simulator::Simulator;
    use std::path::Path;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    // Helper to read all simulation values for a specific site
    fn read_site_simulations(file_path: &Path, site_idx: usize) -> Vec<f64> {
        let file = File::open(file_path).expect("Failed to open CSV");
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        
        // Skip header
        lines.next();
        
        // Skip to site row
        for _ in 0..site_idx {
            lines.next();
        }
        
        let line = lines.next().expect("Site row not found").expect("Read error");
        let parts: Vec<&str> = line.split(',').collect();
        // parts[0] is Site_ID, parts[1..] are Sim_1, Sim_2...
        parts[1..].iter().map(|s| s.parse::<f64>().expect("Parse error")).collect()
    }

    // Helper to calculate median of a vector (Geometric Mean: exp(mean(ln(x))))
    fn calculate_median(vals: Vec<f64>) -> f64 {
        let sum_ln: f64 = vals.iter().map(|v| v.ln()).sum();
        let mean_ln = sum_ln / vals.len() as f64;
        mean_ln.exp()
    }

    // 1. Setup Common Params
    let sites = vec![
        Site::new(1, 120.0, 30.0, 0.0, 0.0, 760.0, None, false),
    ];
    let output_dir = Path::new("output");
    let target_file = output_dir.join("IM_T1.csv");

    // 2. Run Direct Median Simulation (ifmedian=true)
    println!("Running Direct Median Simulation...");
    let eq_source_direct = EQSource::new(
        true, 7.0, 1, 12345, 120.0, 30.0, Some(10.0), Some(20.0), 
        (1.0, 0.0, 0.0), 0.0, false, Some(10.0), Some(0.0), Some(15.0), Region::Global
    );
    let sim_direct = Simulator::new(eq_source_direct, sites.clone(), "CB14".to_string(), 1)
        .with_grid_threshold(1000); // Force direct mode
    sim_direct.run();
    
    let vals_direct = read_site_simulations(&target_file, 0);
    let val_direct = vals_direct[0];

    // 3. Run Random Simulation (ifmedian=false, n_sim=200)
    println!("Running Random Simulation...");
    let eq_source_random = EQSource::new(
        false, 7.0, 200, 12345, 120.0, 30.0, Some(10.0), Some(20.0), 
        (1.0, 0.0, 0.0), 0.0, false, Some(10.0), Some(0.0), Some(15.0), Region::Global
    );
    let sim_random = Simulator::new(eq_source_random, sites.clone(), "CB14".to_string(), 1)
        .with_grid_threshold(1000);
    sim_random.run();

    let vals_random = read_site_simulations(&target_file, 0);
    let val_random_median = calculate_median(vals_random);

    // 4. Compare
    println!("Direct Median: {:.4}, Random Median: {:.4}", val_direct, val_random_median);
    let diff = (val_direct - val_random_median).abs();
    let avg = (val_direct + val_random_median) / 2.0;
    let rel_err = if avg > 1e-9 { diff / avg } else { 0.0 };
    println!("Difference: {:.4}, Rel Error: {:.4}", diff, rel_err);

    assert!(rel_err < 0.02, "Random simulation median should be close to direct median (within 2%)");
}
