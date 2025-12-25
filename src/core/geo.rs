use nalgebra::Vector3;
use proj4rs::proj::Proj;
use proj4rs::transform::transform;

pub type Vec3 = Vector3<f64>;

/// 判断点 p 的垂足是否在四边形 quad 内
/// quad 必须包含 4 个点
pub fn footpoint_inside_quad(p: Vec3, quad: &[Vec3]) -> bool {
    assert_eq!(quad.len(), 4);
    // 矩形两条边
    let p12 = quad[1] - quad[0];
    let p23 = quad[2] - quad[1];
    
    let cross_prod = p12.cross(&p23);
    if cross_prod.norm_squared() < 1e-18 {
        return false;
    }
    
    // 矩形中心
    let c = (quad[0] + quad[1] + quad[2] + quad[3]) * 0.25;
    
    // 法向量
    let mut normal = cross_prod.normalize();
    
    // 确保法向量指向 p 的一侧
    if (c - p).dot(&normal) < 0.0 {
        normal = -normal;
    }
    
    let dist = (c - p).dot(&normal).abs();
    let fp = normal * dist + p;
    
    // 矩形平面内的直线 l = fp - c
    let l = fp - c;
    
    // 相对坐标
    let x = l.dot(&p12.normalize()).abs() / p12.norm();
    let y = l.dot(&p23.normalize()).abs() / p23.norm();
    
    !(x > 0.5 || y > 0.5)
}

/// 计算点 p 到线段 p1-p2 的最小距离
pub fn distance_point_and_line_segment(p: Vec3, p1: Vec3, p2: Vec3) -> f64 {
    let segment_length = (p1 - p2).norm();
    if segment_length < 1e-9 {
        return (p - p1).norm();
    }
    
    // 判断垂足是否在线段外
    if (p - p1).dot(&(p2 - p1)) < 0.0 || (p - p2).dot(&(p1 - p2)) < 0.0 {
        (p - p1).norm().min((p - p2).norm())
    } else {
        (p1 - p).cross(&(p2 - p)).norm() / segment_length
    }
}

// 坐标转换常量
// const R_EARTH: f64 = 6371.393; // 地球半径 km (Deprecated, using proj4rs)

/// 使用等距方位投影 (Azimuthal Equidistant Projection) 将经纬度转换为局部坐标 (km)
/// 
/// 输入:
///  - lon: 经度
///  - lat: 纬度
/// - lon_0: 参考经度 (投影中心)
/// - lat_0: 参考纬度 (投影中心)
///
/// 返回:
/// - (x, y): 局部坐标，单位 km。y 轴指向正北，x 轴指向正东。
pub fn latlon2xy(lon: f64, lat: f64, lon_0: f64, lat_0: f64) -> (f64, f64) {
    // 1. 定义源坐标系 (WGS84 经纬度)
    let wgs84_geo = Proj::from_proj_string(
        "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    ).expect("Failed to create WGS84 projection");

    // 2. 定义目标坐标系 (Azimuthal Equidistant - AEQD)
    // 动态构建投影字符串
    let proj_str = format!(
        "+proj=aeqd +lat_0={} +lon_0={} +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        lat_0, lon_0
    );
    let aeqd_proj = Proj::from_proj_string(&proj_str)
        .expect("Failed to create AEQD projection");

    // 3. 准备输入点 (经度, 纬度, 高度)
    // proj4rs 要求输入必须是弧度
    let mut point = (
        lon.to_radians(), 
        lat.to_radians(), 
        0.0
    );

    // 4. 执行转换
    transform(&wgs84_geo, &aeqd_proj, &mut point)
        .expect("Projection transform failed");

    // 5. 返回结果 (转换为 km)
    (point.0 / 1000.0, point.1 / 1000.0)
}

/// 使用等距方位投影 (Azimuthal Equidistant Projection) 的反变换：
/// 将局部坐标 (km) 转换回经纬度 (度)
///
/// 输入:
/// - x_km, y_km: 局部坐标，单位 km（y 轴指向正北，x 轴指向正东）
/// - lon_0, lat_0: 投影中心经纬度 (度)
///
/// 返回:
/// - (lon, lat): 经纬度 (度)
pub fn xy2latlon(x_km: f64, y_km: f64, lon_0: f64, lat_0: f64) -> (f64, f64) {
    // 源坐标系为 AEQD（单位 m），目标为 WGS84 经纬度（单位弧度）
    let wgs84_geo = Proj::from_proj_string(
        "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    ).expect("Failed to create WGS84 projection");

    let proj_str = format!(
        "+proj=aeqd +lat_0={} +lon_0={} +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        lat_0, lon_0
    );
    let aeqd_proj = Proj::from_proj_string(&proj_str)
        .expect("Failed to create AEQD projection");

    let mut point = (
        x_km * 1000.0,
        y_km * 1000.0,
        0.0
    );

    // AEQD -> WGS84(longlat in radians)
    transform(&aeqd_proj, &wgs84_geo, &mut point)
        .expect("Inverse projection transform failed");

    (point.0.to_degrees(), point.1.to_degrees())
}

/// 计算经度 lon 对应的局部X坐标 (km)，相对于参考经纬度 (lon_0, lat_0)
/// 注意：需要同时提供 lat 以进行准确投影
pub fn get_x(lon: f64, lat: f64, lon_0: f64, lat_0: f64) -> f64 {
    latlon2xy(lon, lat, lon_0, lat_0).0
}


/// 计算纬度 lat 对应的局部Y坐标 (km)，相对于参考经纬度 (lon_0, lat_0)
/// 注意：需要同时提供 lon 以进行准确投影
pub fn get_y(lon: f64, lat: f64, lon_0: f64, lat_0: f64) -> f64 {
    latlon2xy(lon, lat, lon_0, lat_0).1
}

/// 计算 Rrup (Rupture Distance): 场地到断层破裂面的最短距离，单位 km
/// 输入:
/// - site_p: 场地位置，单位 km
/// - rupture_points: 断层破裂面四个顶点位置 (km)，按顺时针或逆时针顺序排列，矩形中心为原点
/// - rupture_normal: 断层面法向量
pub fn calc_rrup(site_p: Vec3, rupture_points: &[Vec3], rupture_normal: Vec3) -> f64 {
    if footpoint_inside_quad(site_p, rupture_points) {
        let normal = rupture_normal.normalize();
        (rupture_points[0] - site_p).dot(&(-normal)).abs()
    } else {
        let d1 = distance_point_and_line_segment(site_p, rupture_points[0], rupture_points[1]);
        let d2 = distance_point_and_line_segment(site_p, rupture_points[1], rupture_points[2]);
        let d3 = distance_point_and_line_segment(site_p, rupture_points[2], rupture_points[3]);
        let d4 = distance_point_and_line_segment(site_p, rupture_points[3], rupture_points[0]);
        d1.min(d2).min(d3).min(d4)
    }
}

/// 计算 Rjb (Joyner-Boore Distance): 场地到断层破裂面在地表投影的最短距离，单位 km
/// 输入:
/// - site_p: 场地位置，单位 km
/// - rupture_points: 断层破裂面四个顶点位置 (km)，按顺时针或逆时针顺序排列
pub fn calc_rjb(site_p: Vec3, rupture_points: &[Vec3]) -> f64 {
    // 投影到 z=0
    let mut rupture_points_proj = rupture_points.to_vec();
    for p in &mut rupture_points_proj {
        p.z = 0.0;
    }
    
    // 场地也投影到 z=0 (虽然传入的 site_p 在调用前通常已经是 z=0，但为了保险起见)
    let mut site_p_proj = site_p;
    site_p_proj.z = 0.0;

    if footpoint_inside_quad(site_p_proj, &rupture_points_proj) {
        0.0
    } else {
        let d1 = distance_point_and_line_segment(site_p_proj, rupture_points_proj[0], rupture_points_proj[1]);
        let d2 = distance_point_and_line_segment(site_p_proj, rupture_points_proj[1], rupture_points_proj[2]);
        let d3 = distance_point_and_line_segment(site_p_proj, rupture_points_proj[2], rupture_points_proj[3]);
        let d4 = distance_point_and_line_segment(site_p_proj, rupture_points_proj[3], rupture_points_proj[0]);
        d1.min(d2).min(d3).min(d4)
    }
}

/// 计算 Rx: 场地到断层迹线（或其延伸线）的水平垂直距离，单位 km
/// 输入:
/// - site_p: 场地位置，单位 km
/// - rupture_points: 断层破裂面四个顶点位置 (km)，按顺时针或逆时针顺序排列
/// - rupture_normal: 断层面法向量
pub fn calc_rx(site_p: Vec3, rupture_points: &[Vec3], rupture_normal: Vec3) -> f64 {
    // 投影到 z=0
    let mut rupture_points_proj = rupture_points.to_vec();
    for p in &mut rupture_points_proj {
        p.z = 0.0;
    }
    let mut site_p_proj = site_p;
    site_p_proj.z = 0.0;
    
    // 垂线距离 (P0-P1 是上边缘)
    let result1 = (rupture_points_proj[0] - site_p_proj).cross(&(rupture_points_proj[1] - site_p_proj)).norm() 
        / (rupture_points_proj[0] - rupture_points_proj[1]).norm();
        
    // 判断 Hanging wall side
    let mut n_xy = rupture_normal;
    n_xy.z = 0.0;
    if n_xy.norm_squared() < 1e-9 {
        panic!("Horizontal fault (normal parallel to Z axis) is not supported.");
    } else {
        n_xy = n_xy.normalize();
    }
    
    let center = (rupture_points_proj[0] + rupture_points_proj[1]) * 0.5;
    if (site_p_proj - center).dot(&n_xy) > 0.0 {
        result1
    } else {
        -result1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-6;

    #[test]
    fn test_get_x_y() {
        // Test origin
        assert!((get_x(100.0, 30.0, 100.0, 30.0) - 0.0).abs() < EPSILON);
        assert!((get_y(100.0, 30.0, 100.0, 30.0) - 0.0).abs() < EPSILON);

        // Test simple offset
        // 1 degree lat is approx 111km
        let y = get_y(100.0, 31.0, 100.0, 30.0);
        // WGS84 distance for 1 degree latitude is approx 110.6 to 111.7 km
        // Old spherical approximation was ~111.2 km
        assert!((y - 111.0).abs() < 1.0); 
    }

    #[test]
    fn test_distance_point_and_line_segment() {
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(10.0, 0.0, 0.0);

        // Point on segment
        let p = Vec3::new(5.0, 0.0, 0.0);
        assert!(distance_point_and_line_segment(p, p1, p2).abs() < EPSILON);

        // Point above segment (projection inside)
        let p = Vec3::new(5.0, 5.0, 0.0);
        assert!((distance_point_and_line_segment(p, p1, p2) - 5.0).abs() < EPSILON);

        // Point outside segment (closer to p1)
        let p = Vec3::new(-3.0, 4.0, 0.0);
        // dist to p1 is 5.0 (3-4-5 triangle)
        assert!((distance_point_and_line_segment(p, p1, p2) - 5.0).abs() < EPSILON);

        // Point outside segment (closer to p2)
        let p = Vec3::new(13.0, 4.0, 0.0);
        // dist to p2 is 5.0
        assert!((distance_point_and_line_segment(p, p1, p2) - 5.0).abs() < EPSILON);
    }

    #[test]
    fn test_footpoint_inside_quad() {
        // Define a 10x10 square on XY plane centered at origin
        // Note: The function expects points in order.
        // p1(5,5), p2(-5,5), p3(-5,-5), p4(5,-5)
        let p1 = Vec3::new(5.0, 5.0, 0.0);
        let p2 = Vec3::new(-5.0, 5.0, 0.0);
        let p3 = Vec3::new(-5.0, -5.0, 0.0);
        let p4 = Vec3::new(5.0, -5.0, 0.0);
        let quad = vec![p1, p2, p3, p4];

        // Point directly above center
        let p = Vec3::new(0.0, 0.0, 10.0);
        assert!(footpoint_inside_quad(p, &quad));

        // Point above edge (inside)
        let p = Vec3::new(4.9, 0.0, 10.0);
        assert!(footpoint_inside_quad(p, &quad));

        // Point outside
        let p = Vec3::new(6.0, 0.0, 10.0);
        assert!(!footpoint_inside_quad(p, &quad));
        
        // Point outside
        let p = Vec3::new(0.0, 6.0, 10.0);
        assert!(!footpoint_inside_quad(p, &quad));
    }

    #[test]
    fn test_xy2latlon_roundtrip() {
        let lon_0 = 100.0;
        let lat_0 = 30.0;
        let lon = 100.2;
        let lat = 30.15;

        let (x, y) = latlon2xy(lon, lat, lon_0, lat_0);
        let (lon_rt, lat_rt) = xy2latlon(x, y, lon_0, lat_0);

        assert!((lon_rt - lon).abs() < 1e-6);
        assert!((lat_rt - lat).abs() < 1e-6);
    }

    #[test]
    fn test_xy2latlon_north_shift() {
        let lon_0 = 100.0;
        let lat_0 = 30.0;
        // 约 111 km 北移应接近 +1 度纬度
        let y_km = 111.0;
        let x_km = 0.0;
        let (lon1, lat1) = xy2latlon(x_km, y_km, lon_0, lat_0);
        assert!((lon1 - lon_0).abs() < 0.01); // 经度变化较小
        assert!((lat1 - (lat_0 + 1.0)).abs() < 0.01); // 近似 1 度（投影误差允许）
    }
}
