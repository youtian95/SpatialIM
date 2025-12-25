use super::site::Site;
use super::geo;
use spade::{DelaunayTriangulation, Point2, Triangulation, PositionInTriangulation};
use std::collections::HashMap;

/// 根据给定场地的经纬度范围，生成正方形网格点（经纬度）。
/// - 网格在局部坐标 (km) 下为正方形，XY 步长一致。
/// - 通过 `grid_spacing_km` 控制网格边长；若为 `None`，默认 0.5 km（500 m）。
/// - 返回经纬度点列表以及网格维度 (nx, ny)。
pub fn generate_grid_points_from_sites(
    sites: &[Site],
    grid_spacing_km: Option<f64>,
) -> (Vec<(f64, f64)>, usize, usize) {
    if sites.is_empty() {
        return (Vec::new(), 0, 0);
    }

    let mut min_lon = sites[0].lon;
    let mut max_lon = sites[0].lon;
    let mut min_lat = sites[0].lat;
    let mut max_lat = sites[0].lat;

    for s in sites.iter() {
        if s.lon < min_lon { min_lon = s.lon; }
        if s.lon > max_lon { max_lon = s.lon; }
        if s.lat < min_lat { min_lat = s.lat; }
        if s.lat > max_lat { max_lat = s.lat; }
    }

    // 投影到 XY (km) 以保证网格为正方形
    let eps = 1e-6;
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    // 默认参考经纬度使用第一个点
    let lon_0 = sites[0].lon;
    let lat_0 = sites[0].lat;

    for s in sites.iter() {
        let (x, y) = geo::latlon2xy(s.lon, s.lat, lon_0, lat_0);
        if x < min_x { min_x = x; }
        if x > max_x { max_x = x; }
        if y < min_y { min_y = y; }
        if y > max_y { max_y = y; }
    }
    if (max_x - min_x).abs() < eps { max_x = min_x + eps; }
    if (max_y - min_y).abs() < eps { max_y = min_y + eps; }

    let range_x = max_x - min_x;
    let range_y = max_y - min_y;

    // 使用指定或默认网格边长（km）
    let d = grid_spacing_km.unwrap_or(0.5).max(eps);

    // 计算网格维度
    // 向上取整，确保网格覆盖整个范围
    let n_intervals_x = (range_x / d).ceil() as usize;
    let n_intervals_y = (range_y / d).ceil() as usize;

    // 计算总覆盖宽度和高度
    let total_w = n_intervals_x as f64 * d;
    let total_h = n_intervals_y as f64 * d;

    // 计算偏移量，使网格居中于范围
    // 这样起始点和终点距离网格两端（min/max）的距离相同
    let offset_x = (total_w - range_x) / 2.0;
    let offset_y = (total_h - range_y) / 2.0;

    // 计算起始坐标
    let start_x = min_x - offset_x;
    let start_y = min_y - offset_y;

    // 网格点数 = 间隔数 + 1
    let nx = n_intervals_x + 1;
    let ny = n_intervals_y + 1;

    // 在 XY 上生成网格点，再用反投影转换为经纬度
    let mut points = Vec::with_capacity(nx * ny);
    for iy in 0..ny {
        let y = start_y + (iy as f64) * d;
        for ix in 0..nx {
            let x = start_x + (ix as f64) * d;
            let (lon, lat) = geo::xy2latlon(x, y, lon_0, lat_0);
            points.push((lon, lat));
        }
    }

    (points, nx, ny)
}

/// 使用 Delaunay 三角插值，将原始场地属性插值到新的网格点。若点落在凸包外则退化为最近邻。周期全部设为 0.0。
/// # 参数：
/// - `sites`: 原始场地
/// - `points`: 新的网格经纬度点 (lon, lat)
/// # 返回：
/// - 插值后的场地列表
pub fn interpolate_sites_to_points(sites: &[Site], points: &[(f64, f64)]) -> Vec<Site> {
    if sites.is_empty() {
        return Vec::new();
    }

    // 投影中心使用第一个场地
    let lon_0 = sites[0].lon;
    let lat_0 = sites[0].lat;

    // 使用 Spade 进行 Delaunay 三角剖分
    // 使用 HashMap 存储 Handle -> Site Index 的映射
    let mut tri: DelaunayTriangulation<Point2<f64>> = DelaunayTriangulation::new();
    let mut site_map = HashMap::new();

    for (i, s) in sites.iter().enumerate() {
        let (x, y) = geo::latlon2xy(s.lon, s.lat, lon_0, lat_0);
        if let Ok(handle) = tri.insert(Point2::new(x, y)) {
            site_map.insert(handle, i);
        }
    }

    let mut result = Vec::with_capacity(points.len());

    for (idx, (lon, lat)) in points.iter().enumerate() {
        let (px, py) = geo::latlon2xy(*lon, *lat, lon_0, lat_0);
        let p = Point2::new(px, py);

        // 定位点在三角剖分中的位置
        let position = tri.locate(p);

        let new_site = match position {
            PositionInTriangulation::OnVertex(v_handle) => {
                let site_idx = site_map[&v_handle];
                let mut s = sites[site_idx].clone();
                s.id = (idx as i32) + 1;
                s.lon = *lon;
                s.lat = *lat;
                s.period1 = 0.0;
                s
            },
            PositionInTriangulation::OnEdge(e_handle) => {
                // 在边上，线性插值
                let edge = tri.directed_edge(e_handle);
                let v1 = edge.from();
                let v2 = edge.to();
                let idx1 = site_map[&v1.fix()];
                let idx2 = site_map[&v2.fix()];
                let p1 = v1.position();
                let p2 = v2.position();
                
                // 计算权重
                let d_total = ((p1.x - p2.x).powi(2) + (p1.y - p2.y).powi(2)).sqrt();
                let d1 = ((p1.x - px).powi(2) + (p1.y - py).powi(2)).sqrt();
                let w2 = if d_total > 1e-9 { d1 / d_total } else { 0.5 };
                let w1 = 1.0 - w2;
                
                interpolate_2_sites(&sites[idx1], &sites[idx2], w1, w2, (idx as i32) + 1, *lon, *lat)
            },
            PositionInTriangulation::OnFace(f_handle) => {
                // 在三角形内，重心插值
                let [v0, v1, v2] = tri.face(f_handle).vertices();
                let idx0 = site_map[&v0.fix()];
                let idx1 = site_map[&v1.fix()];
                let idx2 = site_map[&v2.fix()];
                
                let p0 = v0.position();
                let p1 = v1.position();
                let p2 = v2.position();
                
                // 计算重心坐标
                let (w0, w1, w2) = barycentric_coords(&p0, &p1, &p2, &p);
                
                interpolate_3_sites(&sites[idx0], &sites[idx1], &sites[idx2], w0, w1, w2, (idx as i32) + 1, *lon, *lat)
            },
            PositionInTriangulation::OutsideOfConvexHull(_) | PositionInTriangulation::NoTriangulation => {
                // 凸包外，使用最近邻
                if let Some(v_ref) = tri.nearest_neighbor(p) {
                    let site_idx = site_map[&v_ref.fix()];
                    let mut s = sites[site_idx].clone();
                    s.id = (idx as i32) + 1;
                    s.lon = *lon;
                    s.lat = *lat;
                    s.period1 = 0.0;
                    s
                } else {
                    let mut s = sites[0].clone();
                    s.period1 = 0.0;
                    s
                }
            }
        };
        result.push(new_site);
    }

    result
}

fn barycentric_coords(p0: &Point2<f64>, p1: &Point2<f64>, p2: &Point2<f64>, p: &Point2<f64>) -> (f64, f64, f64) {
    let det = (p1.y - p2.y) * (p0.x - p2.x) + (p2.x - p1.x) * (p0.y - p2.y);
    if det.abs() < 1e-12 {
        return (1.0/3.0, 1.0/3.0, 1.0/3.0);
    }
    let l0 = ((p1.y - p2.y) * (p.x - p2.x) + (p2.x - p1.x) * (p.y - p2.y)) / det;
    let l1 = ((p2.y - p0.y) * (p.x - p2.x) + (p0.x - p2.x) * (p.y - p2.y)) / det;
    let l2 = 1.0 - l0 - l1;
    (l0, l1, l2)
}

fn interpolate_2_sites(s1: &Site, s2: &Site, w1: f64, w2: f64, id: i32, lon: f64, lat: f64) -> Site {
    let interp = |a: f64, b: f64| w1 * a + w2 * b;
    
    let z25 = match (s1.z25, s2.z25) {
        (Some(a), Some(b)) => Some(interp(a, b)),
        _ => None,
    };

    Site {
        id,
        lon,
        lat,
        elevation_km: interp(s1.elevation_km, s2.elevation_km),
        period1: 0.0,
        vs30: interp(s1.vs30, s2.vs30),
        z25,
        r_rup: None,
        r_jb: None,
        r_x: None,
    }
}

fn interpolate_3_sites(s0: &Site, s1: &Site, s2: &Site, w0: f64, w1: f64, w2: f64, id: i32, lon: f64, lat: f64) -> Site {
    let interp = |a: f64, b: f64, c: f64| w0 * a + w1 * b + w2 * c;

    let z25 = match (s0.z25, s1.z25, s2.z25) {
        (Some(a), Some(b), Some(c)) => Some(interp(a, b, c)),
        _ => None,
    };

    Site {
        id,
        lon,
        lat,
        elevation_km: interp(s0.elevation_km, s1.elevation_km, s2.elevation_km),
        period1: 0.0,
        vs30: interp(s0.vs30, s1.vs30, s2.vs30),
        z25,
        r_rup: None,
        r_jb: None,
        r_x: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::site::Site;

    #[test]
    fn test_interpolate_sites_to_points() {
        // 构造 3 个场地形成一个三角形
        // Site 1: (0, 0)
        let s1 = Site::new(1, 0.0, 0.0, 0.0, 1.0, 100.0, Some(1.0), false);
        // Site 2: (1, 0) - approx 111km away
        let s2 = Site::new(2, 1.0, 0.0, 10.0, 2.0, 200.0, Some(2.0), false);
        // Site 3: (0, 1) - approx 111km away
        let s3 = Site::new(3, 0.0, 1.0, 20.0, 3.0, 300.0, Some(3.0), false);

        let sites = vec![s1, s2, s3];

        // 测试点
        let points = vec![
            (0.0, 0.0),   // On Vertex 1
            (0.25, 0.25), // Inside triangle (barycentric)
            (2.0, 2.0),   // Outside convex hull (nearest neighbor)
        ];

        let result = interpolate_sites_to_points(&sites, &points);

        assert_eq!(result.len(), 3);

        // 1. On Vertex 1
        // Should match s1 exactly, except period1 = 0.0
        let r1 = &result[0];
        assert_eq!(r1.id, 1);
        assert!((r1.lon - 0.0).abs() < 1e-6);
        assert!((r1.lat - 0.0).abs() < 1e-6);
        assert!((r1.vs30 - 100.0).abs() < 1e-6);
        assert_eq!(r1.period1, 0.0); // Check period1 is 0.0

        // 2. Inside triangle
        // (0.25, 0.25) is the centroid of (0,0), (0.5,0), (0,0.5) if scaled?
        // Actually for (0,0), (1,0), (0,1), the point (0.25, 0.25) is inside.
        // Barycentric coords?
        // Let's just check that values are between min and max
        let r2 = &result[1];
        assert_eq!(r2.id, 2);
        assert!(r2.vs30 > 100.0 && r2.vs30 < 300.0);
        assert!(r2.elevation_km > 0.0 && r2.elevation_km < 20.0);
        assert_eq!(r2.period1, 0.0); // Check period1 is 0.0

        // 3. Outside convex hull
        // Nearest to (2,2) should be (0,1) or (1,0)?
        // (2,2) is far. (1,0) dist is sqrt(1^2+2^2)=sqrt(5). (0,1) dist is sqrt(2^2+1^2)=sqrt(5).
        // Wait, (1,0) to (2,2) is dx=1, dy=2.
        // (0,1) to (2,2) is dx=2, dy=1.
        // Distances are roughly equal in lat/lon space.
        // But let's check period1 is 0.0 regardless of which one is nearest.
        let r3 = &result[2];
        assert_eq!(r3.id, 3);
        assert_eq!(r3.period1, 0.0); // Check period1 is 0.0
    }
}