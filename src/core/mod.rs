//! # 地震强度场模拟的原理
//! 
//! 统计假设：各场地 IM 服从对数正态分布，记：
//! $$ \ln IM = \mu_{\ln IM} + \tau_{ij}(T) \cdot \eta_{i} + \phi_{ij}(T) \cdot \varepsilon_{ij} $$
//! 其中 $i$ 表示地震事件，$j$ 表示场地，$T$ 表示周期（PGA 可视为 $T=0$）。$\mu_{\ln IM}$ 由经验 GMM（如 NGA‑West2）给出； $\tau_{ij}(T)$、$\phi_{ij}(T)$ 分别为事件间与事件内标准差，均由 GMM 给出，都是跟场地和周期有关系；$\eta_i$ 为事件间标准正态随机变量，在同次地震事件中对所有场地一致；$\varepsilon_{ij}$ 为事件内标准正态随机变量，在不同场地间独立同分布。
//! 
//! 多周期联合模拟流程：
//! 1) 估计均值与方差。使用 NGA‑West2 GMM 估计 $\mu_{\ln IM}$、总标准差 $\sigma(T)$、事件间标准差 $\tau(T)$、事件内标准差$\phi(T)$。
//! 2) 模拟事件间残差 $\delta B$。对于一次地震模拟，所有场地的 $\delta B(T)$相同；不同周期的 $\delta B(T)$ 服从多元正态分布。模拟步骤：
//!     1. 构建跨周期相关矩阵 $\rho_B(0,T_1,T_2)$：
//!         1. 自 Baker & Jayaram (2008) 获取跨周期总相关 $\rho_{\text{total}}(T_1,T_2)$；
//!         1. 自 Loth & Baker (2013) 获取事件内相关 $\rho_W(h,T_1,T_2)$，取 $h=0$ 得 $\rho_W(0,T_1,T_2)$；
//!         1. 对上面公式同一个场地的不同周期IM计算方差，得到下面的公式，可以根据下面这个公式计算未知的 $\rho_B(T_1,T_2)$：
//!             $$ \rho_{\text{total}}(0,T_1,T_2)\ \sigma(T_1)\ \sigma(T_2) = \rho_B(T_1,T_2)\ \tau(T_1)\ \tau(T_2) + \rho_W(0,T_1,T_2)\ \phi(T_1)\ \phi(T_2). $$
//!     1. 基于 $\rho_B(T_1,T_2)$ 构建相关矩阵 $\Sigma_B$，并进行 Cholesky 分解得到下三角矩阵 $\mathbf{L}$；
//!     1. 生成事件间残差向量：
//!         $$ \boldsymbol{\delta B} = \mathbf{L} \boldsymbol{u}, $$
//!         其中 $\boldsymbol{u}$ 为一次独立标准正态随机向量的实现。同一模拟中，所有场地共享同一 $\delta B(T)$ 向量。
//! 3) 模拟空间相关的事件内残差 $\delta W$。事件内残差随周期与空间共同变动，服从多元正态；采用 PCA 等方法高效生成全体场地的相关样本。
//! 4) 组合与输出。对每个场地、每次模拟、每个周期：
//!    $$ S_a(T) = \exp \big( \mu_{\ln IM} + \delta B(T) + \delta W(T) \big). $$
//! 
//! 产出：得到一组随周期的随机 $S_a(T)$ 场（含 PGA），可用于风险评估与结构响应分析，并支持按建筑自振周期插值。
//! 

mod gmpe;
mod b_res_sim;
mod w_res_sim;
pub mod simulator;
pub mod io;
pub mod geo;
pub mod site;
pub mod eq_source;
pub mod utilities;
pub mod grid;


