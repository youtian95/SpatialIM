# 地震动强度空间分布随机场模拟

## 示例
以下为`python/examples`文件夹中一次7级地震下Sa(T=1s)的一次随机模拟分布，矩形框为断层在地面的投影：
![plot_sim1](python/examples/output/plot_sim1.png)
10次模拟的Sa(T=1s)中值分布：
![plot_median](python/examples/output/plot_median.png)

## Python中使用

1. 安装包:
   ```
   pip install spatialim
   ```
1. `python/examples`文件夹包含示例代码`demo.py`，展示了如何使用该包进行地震动强度的空间分布模拟。

## 参考文献

1. K W Campbell, Y Bozorgnia. NGA-West2 Ground Motion Model for the Average Horizontal Components of PGA, PGV, and 5% Damped Linear Acceleration Response Spectra. Earthquake Spectra, 2014, 30(3): 1087-1115.
1. N Jayaram, J W Baker. Correlation model for spatially distributed ground-motion intensities. Earthquake Engineering & Structural Dynamics, 2009, 38(15): 1687-1708.
1. K Goda. Interevent Variability of Spatial Correlation of Peak Ground Motions and Response Spectra. Bulletin of the Seismological Society of America, 2011, 101(5): 2522-2531.
1. M Markhvida, L Ceferino, J W Baker. Modeling spatially correlated spectral accelerations at multiple periods using principal component analysis and geostatistics. Earthquake Engineering & Structural Dynamics, 2018, 47(5): 1107-1123.

## 开发文档

### 代码结构

核心逻辑位于 `src/core/` 目录下：

 - `simulator.rs`: 模拟器，负责调用 GMPE 和残差模拟
 - 功能模块
   - `io.rs`: 处理输入输出文件
   - `geo.rs`: 提供距离计算等功能
   - `site.rs`: 定义场地
   - `eq_source.rs`: 定义震源
   - `utilities.rs`: 提供通用工具函数
 - 核心模块
   - `gmpe`: 提供各种GMPE模型，目前实现了CB14模型
   - `b_res_sim.rs`: 提供模拟地震事件间残差的方法
   - `w_res_sim.rs`: 提供模拟地震事件内残差的方法

### 分发为Python包

1. 安装 `maturin`:
   ```
   pip install maturin
   ```
1. 修改`Cargo.toml`，添加以下内容:
   ```toml
   [lib]
   name = "_spatialim"
   crate-type = ["cdylib", "rlib"]

   [dependencies]
   pyo3 = { version = "0.27.1", features = ["extension-module"] }
   ```
1. 添加 `pyproject.toml` 文件:
   ```toml
   [build-system]
   requires = ["maturin>=1.0,<2.0"]
   build-backend = "maturin"

   [tool.maturin]
   python-source = "python"
   module-name = "spatialim._spatialim"
   exclude = ["**/*.pyd", "**/*.so", "**/*.dylib"]
   ```
   其他内容与常规 `pyproject.toml` 文件相同。
1. 创建 `python` 目录，结构如下:
   ```
   python/
   ├── spatialim/
   │   ├── __init__.py
   │   ├── _spatialim.pyi # 类型提示文件
   │   └── 其他模块文件.py
   └── setup.py
   ```
1. 开发调试，将包安装到当前Python环境:
   ```
   maturin develop --release
   ```
1. 构建包:
   ```
   maturin build --release
   ```


