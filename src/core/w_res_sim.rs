//! 事件内残差模拟


use nalgebra::{DMatrix, DVector, Cholesky, SymmetricEigen};
use rand::prelude::*;
use rand_distr::{Normal, Distribution};

/// 变差函数模型结构体 / Variogram model struct
#[derive(Debug, Clone)]
struct ModelVario {
    cn: f64,
    c1: f64,
    a1: f64,
    c2: f64,
    a2: f64,
    type_: String, // "nug" or "iso nest"
}

impl ModelVario {
    /// 创建 Nugget 模型
    fn new_nug(cn: f64) -> Self {
        Self {
            cn,
            c1: 0.0,
            a1: 0.0,
            c2: 0.0,
            a2: 0.0,
            type_: "nug".to_string(),
        }
    }

    /// 创建 Isotropic Nested 模型
    fn new_iso_nest(cn: f64, c1: f64, a1: f64, c2: f64, a2: f64) -> Self {
        Self {
            cn,
            c1,
            a1,
            c2,
            a2,
            type_: "iso nest".to_string(),
        }
    }
}

/// 事件内残差模拟器
/// 
/// 基于 PCA 和地统计学方法模拟空间相关的地震动强度残差。
/// 参考文献: Markhvida M, Ceferino L, Baker JW (2018)
pub struct WResSimulator {
    /// 主成分数量
    n_pcs: usize,
    /// 模拟周期列表
    t_periods: DVector<f64>,
    /// PCA 系数矩阵
    pca_coefs: DMatrix<f64>,
    /// 方差缩放因子
    variance_scale_factor: DVector<f64>,
    /// 变差函数模型列表
    model_vario: Vec<ModelVario>,
    /// 分解后的协方差矩阵 (L 矩阵)
    l_matrices: Vec<DMatrix<f64>>, 
    /// 随机数生成器
    rng: StdRng,
}

impl WResSimulator {
    /// 创建新的模拟器实例
    /// 
    /// 初始化 PCA 系数、变差函数参数等固定数据。
    pub fn new(n_pcs: Option<usize>, seed: u64) -> Self {
        let n_pcs = n_pcs.unwrap_or(6);
        let t_periods = DVector::from_vec(vec![
            0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0
        ]);

        #[rustfmt::skip]
        let pca_coefs_data = vec![
            0.270963956240108, -0.139418157111539, 0.0690420060759825, -0.106094866059058, -0.0922880747536406, -0.11348997612886, -0.188935371416588, 0.153956801827048, -0.160082932478587, -0.0485878661500656, 0.106169114128661, 0.0545367125260001, -0.0842347288819985, 0.00206507178166717, 0.233666515544382, -0.0444106080542835, -0.298766213268524, -0.527588859528322, -0.580349072958488,
            0.270185457409379, -0.141734438837191, 0.0770156687411159, -0.11639353409999, -0.103464378298138, -0.124082463392091, -0.199840301009291, 0.155452551301531, -0.157024101304166, -0.0511781532022337, 0.102685985245256, 0.053409178198703, -0.0785807732589693, 0.00538047460862058, 0.220317828901249, -0.039452593127257, -0.257172668625987, -0.15099488917774, 0.781868928002583,
            0.266716484131893, -0.150918021372557, 0.101241750461614, -0.144620230365225, -0.128327845057672, -0.150413273486784, -0.21751911461417, 0.15453312842206, -0.14455513348748, -0.0493784913365888, 0.0865380808702005, 0.0369034473261574, -0.0551975904000634, 0.00787249481774084, 0.149651509850332, -0.0232087969820605, -0.0284597879541604, 0.808901723567414, -0.226437334732122,
            0.251688240452975, -0.184642998841202, 0.178879968100826, -0.221328310597117, -0.17555752575313, -0.176668887066346, -0.1886513513666, 0.0424749058524636, -0.0455090101766787, -0.0291885726697761, -0.0315764520880226, -0.0605411454951765, 0.0935508235649957, 0.0224423933374744, -0.299181352126819, 0.0599168858817312, 0.754350402692246, -0.206472972801382, 0.023109344507298,
            0.236434540660266, -0.218922079202474, 0.23725418394185, -0.234559034122274, -0.133267087790244, -0.0431828093963935, 0.119447151365187, -0.272310554909179, 0.238192698302355, 0.100676333288767, -0.263315034080161, -0.121207177353811, 0.202769694434657, 0.00661936093850581, -0.493306767118467, 0.116246172684757, -0.475923822815555, 0.0367733765077595, -0.00643955732220411,
            0.232994643271854, -0.228087987254185, 0.230554572919478, -0.16044302411133, 0.0400564218872398, 0.181657484726543, 0.427112684239685, -0.324579279763868, 0.263780433255153, 0.142634796459082, -0.0813780347714383, 0.0465305509298986, -0.151801546733413, -0.0833183140197916, 0.534198178434453, -0.184596749983802, 0.210357916588565, -0.0028422556361022, 0.00331108666921881,
            0.238919759244457, -0.211905954003063, 0.132646222385762, 0.0820453502922968, 0.32794697293789, 0.393273105011282, 0.325836316409324, 0.162029624546621, -0.182164846060428, -0.138319895254022, 0.470111475270552, 0.177876087511963, -0.1112565000941, 0.0883177907230135, -0.291143253148395, 0.262494244929844, -0.00152291509197, 0.0154770015927024, 0.00150080927064067,
            0.247247200513419, -0.174053609784519, -0.00825743327819653, 0.277382297057791, 0.403271333843203, 0.220437620135438, -0.0837312940200531, 0.224796020087722, -0.171941472633753, -0.0292747130158189, -0.38152412132245, -0.237244495814365, 0.356271008578838, -0.0850494996785505, -0.0125434811410974, -0.442559354727306, 0.0151125187962444, 0.0110964548458275, 0.00194187520085148,
            0.253677096925723, -0.122375885147353, -0.148595585559558, 0.365271223153003, 0.253186443805678, -0.0612442510569459, -0.283389735157485, -0.0811437931924326, 0.212210668293107, 0.14336242671851, -0.275727618896954, -0.0411307164923357, -0.202014525375018, 0.022845903851548, 0.155257213551213, 0.632145331896681, 0.0455130187869974, 0.000690596761983894, 0.000469634928915885,
            0.254921191707501, -0.071319446353373, -0.237030888421264, 0.359073100317565, 0.0401080106567685, -0.248766601944991, -0.14185903779119, -0.286692239119504, 0.300971237941779, 0.0579993716202019, 0.328411991836789, 0.208361703130758, -0.194768507716927, 0.0324295008946133, -0.258822160901248, -0.477244327079575, 0.00191094743796522, 0.00662185258908696, 0.000188266912551749,
            0.252458254214951, 0.0125091293591534, -0.32712108086482, 0.226053196913636, -0.26129762020473, -0.216236975254179, 0.344080559560279, -0.121230620609225, -0.0602714322805403, -0.219189670580381, 0.211470671417798, -0.128634841134357, 0.576234521196323, -0.0550760430415841, 0.197333043823494, 0.20576645538163, 0.023662144998231, 0.00370370109158968, 0.000185197342559174,
            0.24594424065373, 0.0799604139803053, -0.358449872812475, 0.064099810496083, -0.341792253779399, 0.022496754565775, 0.388717982330275, 0.177122203733985, -0.255990758059885, -0.00644303562267801, -0.37539012296473, -0.0762061466572891, -0.502002035655987, 0.0183662355165662, -0.176351451500561, -0.0686345486566678, 0.0154233464320777, 0.00541822807684994, 0.00127868850390311,
            0.225758567264058, 0.191381035473567, -0.335176303224685, -0.216152633771253, -0.165178954634192, 0.423011571619087, -0.144255461671014, 0.187567292296388, 0.149360081833295, 0.530105300516035, 0.041796232067731, 0.326784764375549, 0.27460983615727, 0.0558724755215897, 0.00442737635788848, 0.0111352503084675, 0.0244045338837114, 0.000682889258092345, 0.00019866537927149,
            0.211097168683674, 0.259405649803992, -0.243643584807687, -0.325745718814885, 0.0763484285808871, 0.330279571319951, -0.220001696329049, -0.117395307490011, 0.271296590718256, -0.438774328341411, 0.148165305851337, -0.48454567943612, -0.143691433970383, -0.0389147466608019, 0.00893381613114769, -0.0205549231060488, -0.00611584426449596, 0.00380246484862132, -0.000389229000503269,
            0.188387436860863, 0.329799740851314, -0.0946692611641819, -0.273646378894757, 0.356651426359913, -0.153161129681719, -0.000682050969765767, -0.329897663448327, -0.267361749795915, -0.2793821559276, -0.26374050060882, 0.528987613257165, 0.0703281456138742, -0.0835258438506085, -0.0254953444060441, 0.0271139318663258, 0.0126194670232453, 0.00385357399236438, -0.000809502038431328,
            0.176395533106338, 0.357332294420881, 0.0554387508851632, -0.155161957798969, 0.354513035089274, -0.343041979311133, 0.161895880026785, -0.0275355960813715, -0.20785940912098, 0.506833313170992, 0.205450556183935, -0.413916086467618, -0.0407497250821277, 0.168472840930613, -0.00235465742383666, -0.00588317616785561, -0.00334222371106789, 0.00197868686362638, 0.00170392027445259,
            0.165469018345669, 0.360040619637271, 0.260392170105291, 0.0669803672081155, 0.0572701975680616, -0.220913199055637, 0.181062550106522, 0.519913777311151, 0.462086625105155, -0.104489884552583, -0.0219632037127792, 0.119011798891668, -0.00429988165275597, -0.417545575417374, -0.0401081045672218, 0.020060720387271, -0.00526025219733044, -0.00580167364426092, 0.00045823003408634,
            0.159580892256856, 0.347927159738324, 0.348346295207843, 0.239394140984691, -0.157211928082521, 0.0928779108728929, -0.00501602103568442, 0.0169759687957175, 0.109978222336614, -0.182603153808833, -0.121233489669211, 0.0711306930722663, 0.0620582733857731, 0.750450107378881, 0.0785420660886401, -0.0521102089376762, 0.0075193664688849, -0.00188268149870413, -0.00185190445021904,
            0.1488329207305, 0.33284769539218, 0.36509805168606, 0.3312276948947, -0.281334582456685, 0.283334016381439, -0.182848344502475, -0.325817997238157, -0.310946153889112, 0.128556113560638, 0.0837173663270596, -0.0703828760736262, -0.046192412575096, -0.441499963638099, -0.0405906261145499, 0.033716161820098, 0.00116291918248088, 0.0040160503301927, 0.000505570344725481
        ];
        let pca_coefs = DMatrix::from_row_slice(19, 19, &pca_coefs_data);

        let variance_scale_factor = DVector::from_vec(vec![
            0.639841744964052, 0.846277138923352, 0.904533062923065, 0.933402823301903, 0.950154585181316, 0.960461560434709, 0.967972143777431, 0.973878205021816, 0.979297028245927, 0.983345801486615, 0.986639491844098, 0.989681397350252, 0.992367578022788, 0.994790440764041, 0.996929078974732, 0.998759738778621, 0.999764879971799, 0.999969807143844, 1.0
        ]);

        let mut model_vario = Vec::with_capacity(19);
        model_vario.push(ModelVario::new_iso_nest(2.500000000000001, 4.520000000000002, 15.0, 6.780000000000003, 250.0));
        model_vario.push(ModelVario::new_iso_nest(0.500000000000000, 1.400000000000000, 10.0, 2.600000000000001, 160.0));
        model_vario.push(ModelVario::new_iso_nest(0.150000000000000, 0.420000000000000, 15.0, 0.630000000000000, 160.0));
        model_vario.push(ModelVario::new_iso_nest(0.150000000000000, 0.225000000000000, 10.0, 0.225000000000000, 120.0));
        model_vario.push(ModelVario::new_nug(0.314321867136085));
        model_vario.push(ModelVario::new_nug(0.190749535511365));
        model_vario.push(ModelVario::new_nug(0.137846758971697));
        model_vario.push(ModelVario::new_nug(0.111283843493347));
        model_vario.push(ModelVario::new_nug(0.096499281204429));
        model_vario.push(ModelVario::new_nug(0.071736796680044));
        model_vario.push(ModelVario::new_nug(0.064816215163269));
        model_vario.push(ModelVario::new_nug(0.054076635653157));
        model_vario.push(ModelVario::new_nug(0.051188751201166));
        model_vario.push(ModelVario::new_nug(0.043316419835114));
        model_vario.push(ModelVario::new_nug(0.041398046055658));
        model_vario.push(ModelVario::new_nug(0.034663671578289));
        model_vario.push(ModelVario::new_nug(0.018796994070527));
        model_vario.push(ModelVario::new_nug(0.002856941187582));
        model_vario.push(ModelVario::new_nug(3.606453961200486e-04));

        Self {
            n_pcs,
            t_periods,
            pca_coefs,
            variance_scale_factor,
            model_vario,
            l_matrices: Vec::new(),
            rng: StdRng::seed_from_u64(seed),
        }
    }

    /// 注册场地并计算空间相关性
    /// 
    /// 计算场地间的距离矩阵，并根据距离矩阵和变差函数模型构建协方差矩阵。
    /// 对协方差矩阵进行 Cholesky 分解（或特征值分解作为回退），用于后续的随机场生成。
    pub fn register_sites(&mut self, x: &DVector<f64>, y: &DVector<f64>) {
        let _n_locs = x.len();
        let distance_matrix = Self::get_distance_matrix(x, y);

        // Scale variance if less than 19 principal components are used
        if self.n_pcs < 19 {
            for i in 0..self.model_vario.len() {
                if self.model_vario[i].type_ == "nug" {
                     self.model_vario[i].cn /= self.variance_scale_factor[self.n_pcs - 1];
                } else {
                     let scale = self.variance_scale_factor[self.n_pcs - 1]; // Index nPCs - 1
                     self.model_vario[i].cn /= scale;
                     self.model_vario[i].c1 /= scale;
                     self.model_vario[i].c2 /= scale;
                }
            }
        }

        self.l_matrices = Vec::with_capacity(self.n_pcs);

        for i in 0..self.n_pcs {
            let cov_matrix = Self::get_covariance(&distance_matrix, &self.model_vario[i]);
            
            // Cholesky decomposition
            let l_matrix = match Cholesky::new(cov_matrix.clone()) {
                Some(cholesky) => cholesky.l(),
                None => {
                    // Fallback to Eigen decomposition
                    let eigen = SymmetricEigen::new(cov_matrix);
                    let eigenvectors = eigen.eigenvectors;
                    let eigenvalues = eigen.eigenvalues;
                    // normTransform = eigenvectors * eigenvalues.sqrt().as_diagonal()
                    let sqrt_eigenvalues = eigenvalues.map(|v| v.sqrt());
                    eigenvectors * DMatrix::from_diagonal(&sqrt_eigenvalues)
                }
            };
            self.l_matrices.push(l_matrix);
        }
    }

    /// 模拟残差
    /// 
    /// 生成指定数量 (n_sims) 的空间相关残差场。
    /// 
    /// # 参数
    /// * `t_sim` - 需要模拟的周期列表
    /// * `n_sims` - 模拟次数
    /// * `phi_matrix` - 可选的事件内标准差矩阵 [n_locs x n_periods]，用于缩放残差。如果没有提供，那么输出的残差的标准差为 1。
    /// 
    /// # 返回
    ///  - `Vec<DMatrix<f64>>`，其中每个矩阵代表一次模拟结果。
    ///     矩阵维度为 [n_locs x n_periods]。
    pub fn simulate_residuals(
        &mut self, 
        t_sim: &DVector<f64>, 
        n_sims: usize,
        phi_matrix: Option<&DMatrix<f64>>
    ) -> Vec<DMatrix<f64>> {
        let n_locs = self.l_matrices[0].nrows();
        
        // Simulate each of the PC's
        // sim_PCA: Vector of Matrices [nLocs x nsims]
        let mut sim_pca = Vec::with_capacity(self.n_pcs);
        let normal = Normal::new(0.0, 1.0).unwrap();

        for i_pc in 0..self.n_pcs {
            // L * randN
            // We can do this efficiently by generating a large random matrix
            // Random matrix R [nLocs x nsims]
            // Samples = L * R
            
            let mut r = DMatrix::zeros(n_locs, n_sims);
            for i in 0..r.len() {
                r[i] = normal.sample(&mut self.rng);
            }
            
            let samples = &self.l_matrices[i_pc] * r;
            sim_pca.push(samples);
        }

        // Transform simulated PC's to spectral acceleration residuals
        let mut sim_results = Vec::with_capacity(n_sims);
        for _ in 0..n_sims {
            sim_results.push(DMatrix::zeros(n_locs, t_sim.len()));
        }

        for (i, &t_val) in t_sim.iter().enumerate() {
            // Find if t_val is in self.t_periods
            let index_opt = self.t_periods.iter().position(|&x| (x - t_val).abs() < 1e-6);

            for j in 0..n_sims {
                // Construct temp_sim_pca for this simulation j: [nLocs x nPCs]
                let mut temp_sim_pca = DMatrix::zeros(n_locs, self.n_pcs);
                for col in 0..self.n_pcs {
                    let col_vec = sim_pca[col].column(j);
                    temp_sim_pca.set_column(col, &col_vec);
                }

                let result_col = if let Some(idx) = index_opt {
                    // Exact match
                    let coefs = self.pca_coefs.row(idx).columns(0, self.n_pcs).transpose(); // [nPCs x 1]
                    // PCA_T_X: transformed * coef.transpose() + mu
                    // Here coef is row vector in C++, passed as matrix.
                    // C++: PCA_T_X(temp_sim_PCA, temp, Zero)
                    // temp is 1xNPCs.
                    // PCA_T_X: original = transformed * coef.transpose()
                    // [nLocs x nPCs] * [nPCs x 1] = [nLocs x 1]
                    temp_sim_pca * coefs
                } else {
                    // Interpolation
                    let mut extra_pca_coefs = DVector::zeros(self.n_pcs);
                    for k in 0..self.n_pcs {
                        let pca_coefs_t = self.pca_coefs.column(k).into_owned();
                        extra_pca_coefs[k] = Self::interp1(&self.t_periods, &pca_coefs_t, t_val);
                    }
                    // extra_pca_coefs is [nPCs x 1] (effectively, though DVector is col vector)
                    temp_sim_pca * extra_pca_coefs
                };

                sim_results[j].set_column(i, &result_col);
            }
        }

        // Apply phi scaling if provided
        if let Some(phi) = phi_matrix {
            for res in &mut sim_results {
                // res is [n_locs x n_periods]
                // phi is [n_locs x n_periods]
                // Element-wise multiplication
                res.component_mul_assign(phi);
            }
        }

        sim_results
    }

    /// 计算距离矩阵
    fn get_distance_matrix(x: &DVector<f64>, y: &DVector<f64>) -> DMatrix<f64> {
        let n = x.len();
        let mut dist = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let d = ((x[i] - x[j]).powi(2) + (y[i] - y[j]).powi(2)).sqrt();
                dist[(i, j)] = d;
            }
        }
        dist
    }

    /// 计算协方差矩阵
    fn get_covariance(dist: &DMatrix<f64>, model: &ModelVario) -> DMatrix<f64> {
        if model.type_ == "iso nest" {
            Self::get_iso_nested_cov(model, dist)
        } else {
            Self::get_nug_cov(model, dist)
        }
    }

    fn get_iso_nested_cov(model: &ModelVario, dist: &DMatrix<f64>) -> DMatrix<f64> {
        let var = model.cn + model.c1 + model.c2;
        let mut cov = DMatrix::zeros(dist.nrows(), dist.ncols());
        
        for i in 0..dist.len() {
            let h = dist[i];
            if h == 0.0 {
                cov[i] = var;
            } else {
                let term1 = model.c1 * (1.0 - (-3.0 * h / model.a1).exp());
                let term2 = model.c2 * (1.0 - (-3.0 * h / model.a2).exp());
                cov[i] = var - (model.cn + term1 + term2);
            }
        }
        cov
    }

    fn get_nug_cov(model: &ModelVario, dist: &DMatrix<f64>) -> DMatrix<f64> {
        let mut cov = DMatrix::zeros(dist.nrows(), dist.ncols());
        for i in 0..dist.len() {
            if dist[i] == 0.0 {
                cov[i] = model.cn;
            } else {
                cov[i] = 0.0;
            }
        }
        cov
    }

    /// 线性插值函数
    fn interp1(x: &DVector<f64>, y: &DVector<f64>, vx: f64) -> f64 {
        if vx < x[0] {
            y[0] + (y[1] - y[0]) / (x[1] - x[0]) * (vx - x[0])
        } else if vx > x[x.len() - 1] {
            let n = x.len();
            y[n - 1] + (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]) * (vx - x[n - 1])
        } else {
            let mut i = 0;
            while i < x.len() - 1 {
                if x[i] <= vx && x[i + 1] >= vx {
                    break;
                }
                i += 1;
            }
            y[i] + (y[i + 1] - y[i]) / (x[i + 1] - x[i]) * (vx - x[i])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    /// 测试生成相关性曲线，需要自己运行后检查结果文件跟文献对比
    #[test]
    fn test_generate_correlation_curves() {
        let periods = vec![0.01, 0.1, 1.0];
        let n_pcs = 19; // Use full model to match literature/MATLAB
        let n_sims = 5000;
        
        // Define sites on a line
        // 1km spacing -> 0 to 200km
        let n_sites = 201;
        let spacing = 1.0;
        let mut x = DVector::zeros(n_sites);
        let mut y = DVector::zeros(n_sites);
        for i in 0..n_sites {
            x[i] = i as f64 * spacing;
            y[i] = 0.0;
        }

        let t_sim = DVector::from_vec(periods.clone());

        let mut sim = WResSimulator::new(Some(n_pcs), 12345);
        sim.register_sites(&x, &y);
        
        // Run simulation
        // results: Vec<DMatrix> of size n_sims. Each matrix is [n_sites x n_periods_sim]
        let results = sim.simulate_residuals(&t_sim, n_sims, None);
        
        // Pre-organize data: [n_periods x n_sites x n_sims]
        let mut period_data = Vec::new();
        for p_idx in 0..periods.len() {
            let mut site_data = Vec::new();
            for s_idx in 0..n_sites {
                let mut vals = Vec::with_capacity(n_sims);
                for k in 0..n_sims {
                    vals.push(results[k][(s_idx, p_idx)]);
                }
                site_data.push(vals);
            }
            period_data.push(site_data);
        }

        // Pre-process data for faster correlation calculation
        // centered_data: [n_periods][n_sites][n_sims]
        // norms: [n_periods][n_sites]
        let mut centered_data = Vec::new();
        let mut norms = Vec::new();

        for p_idx in 0..periods.len() {
            let mut site_centered = Vec::new();
            let mut site_norms = Vec::new();
            
            for s_idx in 0..n_sites {
                let vals = &period_data[p_idx][s_idx];
                let mean = vals.iter().sum::<f64>() / n_sims as f64;
                let mut centered = Vec::with_capacity(n_sims);
                let mut sum_sq = 0.0;
                for &v in vals {
                    let c = v - mean;
                    centered.push(c);
                    sum_sq += c * c;
                }
                site_centered.push(centered);
                site_norms.push(sum_sq.sqrt());
            }
            
            centered_data.push(site_centered);
            norms.push(site_norms);
        }

        let mut file_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        file_path.push("tests/test_result_within_event_correlation.csv");
        let mut file = File::create(file_path).unwrap();
        
        let mut header = "Distance".to_string();
        for &t1 in &periods {
            for &t2 in &periods {
                header.push_str(&format!(",T{}_vs_T{}", t1, t2));
            }
        }
        writeln!(file, "{}", header).unwrap();

        for lag in 0..n_sites {
            let dist = lag as f64 * spacing;
            let mut row = format!("{:.4}", dist);

            for idx1 in 0..periods.len() {
                for idx2 in 0..periods.len() {
                    let mut sum_rho = 0.0;
                    let mut count = 0;

                    // Average over all pairs separated by `lag`
                    for i in 0..(n_sites - lag) {
                        let j = i + lag;
                        
                        let v1 = &centered_data[idx1][i];
                        let v2 = &centered_data[idx2][j];
                        let norm1 = norms[idx1][i];
                        let norm2 = norms[idx2][j];

                        let mut num = 0.0;
                        for k in 0..n_sims {
                            num += v1[k] * v2[k];
                        }
                        
                        if norm1 > 0.0 && norm2 > 0.0 {
                            sum_rho += num / (norm1 * norm2);
                            count += 1;
                        }
                    }
                    
                    let avg_rho = if count > 0 { sum_rho / count as f64 } else { 0.0 };
                    row.push_str(&format!(",{:.6}", avg_rho));
                }
            }
            writeln!(file, "{}", row).unwrap();
        }
    }

    /// 测试生成半变异函数曲线 (Semivariogram)
    /// Gamma(h) = 0.5 * E[(Z(x) - Z(x+h))^2]
    #[test]
    fn test_generate_semivariogram() {
        let periods = vec![0.01, 0.1, 1.0];
        let n_pcs = 19;
        let n_sims = 5000;
        
        // Define sites on a line
        let n_sites = 201;
        let spacing = 1.0;
        let mut x = DVector::zeros(n_sites);
        let mut y = DVector::zeros(n_sites);
        for i in 0..n_sites {
            x[i] = i as f64 * spacing;
            y[i] = 0.0;
        }

        let t_sim = DVector::from_vec(periods.clone());

        let mut sim = WResSimulator::new(Some(n_pcs), 12345);
        sim.register_sites(&x, &y);
        
        let results = sim.simulate_residuals(&t_sim, n_sims, None);
        
        // Organize data: [n_periods][n_sites][n_sims]
        let mut period_data = Vec::new();
        for p_idx in 0..periods.len() {
            let mut site_data = Vec::new();
            for s_idx in 0..n_sites {
                let mut vals = Vec::with_capacity(n_sims);
                for k in 0..n_sims {
                    vals.push(results[k][(s_idx, p_idx)]);
                }
                site_data.push(vals);
            }
            period_data.push(site_data);
        }

        let mut file_path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        file_path.push("tests/test_result_within_event_semivariogram.csv");
        let mut file = File::create(file_path).unwrap();
        
        let mut header = "Distance".to_string();
        for &t in &periods {
            header.push_str(&format!(",T{}", t));
        }
        writeln!(file, "{}", header).unwrap();

        for lag in 0..n_sites {
            let dist = lag as f64 * spacing;
            let mut row = format!("{:.4}", dist);

            for p_idx in 0..periods.len() {
                let mut sum_sq_diff = 0.0;
                let mut count = 0;

                for i in 0..(n_sites - lag) {
                    let j = i + lag;
                    
                    let v1 = &period_data[p_idx][i];
                    let v2 = &period_data[p_idx][j];

                    for k in 0..n_sims {
                        let diff = v1[k] - v2[k];
                        sum_sq_diff += diff * diff;
                    }
                    count += n_sims;
                }
                
                let gamma = if count > 0 { 
                    0.5 * sum_sq_diff / count as f64 
                } else { 
                    0.0 
                };
                row.push_str(&format!(",{:.6}", gamma));
            }
            writeln!(file, "{}", row).unwrap();
        }
    }
}
