# 地震动强度空间分布随机场模拟

## 示例
以下为一次7级地震下PGA的分布，矩形框为断层在地面的投影。例子文件：`Examples\Example 2 - python import\example_plot_IM.py`。
![PGA_contour](./Figures/PGA_contour_M7.png)

## Python中使用

`python/examples`文件夹包含示例

### 使用方法

1. 安装:
   ```
   pip install spatialim
   ```
1. 导入:
   ```python
   import spatialim
   ```

## 代码结构

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

## 参考文献

[1] K W Campbell, Y Bozorgnia. NGA-West2 Ground Motion Model for the Average Horizontal Components of PGA, PGV, and 5% Damped Linear Acceleration Response Spectra. Earthquake Spectra, 2014, 30(3): 1087-1115.

[2] N Jayaram, J W Baker. Correlation model for spatially distributed ground-motion intensities. Earthquake Engineering & Structural Dynamics, 2009, 38(15): 1687-1708.

[3] K Goda. Interevent Variability of Spatial Correlation of Peak Ground Motions and Response Spectra. Bulletin of the Seismological Society of America, 2011, 101(5): 2522-2531.

[4] M Markhvida, L Ceferino, J W Baker. Modeling spatially correlated spectral accelerations at multiple periods using principal component analysis and geostatistics. Earthquake Engineering & Structural Dynamics, 2018, 47(5): 1107-1123.
