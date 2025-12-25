
/// 线性插值
/// y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
pub fn linear_interp(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    if (x1 - x0).abs() < 1e-9 {
        return y0;
    }
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

/// 通用插值（不要求 `xs` 递增），按以下规则：
/// - 找到不递增数组 `xs` 中在目标 `x` 左右最接近的两个点（分别为 `x_left<=x` 与 `x_right>=x`，距离最近）
/// - 在线性插值，并返回插值结果
/// - 若 `x` 超过上下限，则直接返回上下限对应的值（不外推）
/// - 若恰好命中某一点（左、右相同），返回该点的值
pub fn interp_clamped_unsorted(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    if xs.is_empty() || ys.is_empty() || xs.len() != ys.len() {
        return 0.0;
    }

    // 找到上下限（最小与最大 x）
    let mut min_idx = 0usize;
    let mut max_idx = 0usize;
    for i in 1..xs.len() {
        if xs[i] < xs[min_idx] {
            min_idx = i;
        }
        if xs[i] > xs[max_idx] {
            max_idx = i;
        }
    }

    // 边界夹持：超下限或超上限直接返回对应值
    if x <= xs[min_idx] {
        return ys[min_idx];
    }
    if x >= xs[max_idx] {
        return ys[max_idx];
    }

    // 找最近的左点（x_left <= x 且 x_left 最大）与右点（x_right >= x 且 x_right 最小）
    let mut left_idx: Option<usize> = None;
    let mut right_idx: Option<usize> = None;
    let mut best_left_val = f64::NEG_INFINITY;
    let mut best_right_val = f64::INFINITY;

    for i in 0..xs.len() {
        let xi = xs[i];
        if xi <= x && xi > best_left_val {
            best_left_val = xi;
            left_idx = Some(i);
        }
        if xi >= x && xi < best_right_val {
            best_right_val = xi;
            right_idx = Some(i);
        }
    }

    // 理论上在边界排除后应当能找到左右点
    match (left_idx, right_idx) {
        (Some(li), Some(ri)) => {
            if xs[li] == xs[ri] {
                return ys[li];
            }
            return linear_interp(xs[li], ys[li], xs[ri], ys[ri], x);
        }
        (Some(li), None) => ys[li], // 兜底：若没右点，返回左点值
        (None, Some(ri)) => ys[ri], // 兜底：若没左点，返回右点值
        _ => ys[min_idx],           // 进一步兜底
    }
}
