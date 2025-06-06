# 地震动强度空间分布随机场模拟

## 示例
以下为一次7级地震下PGA的分布，矩形框为断层在地面的投影。例子文件：`Examples\Example 2 - python import\example_plot_IM.py`。
![PGA_contour](./Figures/PGA_contour_M7.png)

## 使用

### 示例1 - 可执行文件

Examples/Example 1文件夹中为例子，放入两个文件，名字为EQSource.txt，SiteFile.txt，依次震源信息和场地信息的文件名，然后直接运行IMSim.exe程序，即可进行模拟。

### 输入

1. **EQSource.txt** 每行依次为（每行内不同参数用空格分开）：
    - ifmedian - 0/1，是否输出中位值
    - M - 震级
    - N_sim - 次数
    - seed - int, 随机数种子
    - lon_0，lat_0 - 震中经纬度，°
    - W - 断层破裂面矩形的宽度，km，未知时可输入 999
    - length - 断裂面矩形的长度，km
    - RuptureNormal_x, RuptureNormal_y, RuptureNormal_z - 断裂面朝上的法线方向（向东为x,向北为y,向上为z）
    - lambda - 走滑角（°）- 上盘在破裂面内测量的滑移平均角度，与strike方向相同为0度，逆时针为正值
    - Fhw - 是否考虑上盘效应，0/1
    - Zhyp - 从海平面测量的震源深度，km， unknown, 未知时可输入 999
    - region - 研究的区域
     = 0 全球 (包括台湾)
     = 1 加州
     = 3 中国或者土耳其
     = 4 意大利
    - nPCs - IM相关性PCA方法模拟考虑的主成分阶数，推荐大于等于5
1. **SiteFile.txt** 每行为一个场地的数据，每一行空格分开依次为
    - ID - 场地点的编号
    - lon - 经度
    - lat - 纬度
    - elevation_km - 高程，km
    - period1 - 基本周期
    - Vs30_mpers - 剪切波速
    - Z25_km - 到2.5km/s剪切波速水平面的深度，km，（如果在加州或者日本， Z25_km未知, 可以输入999）

### 输出

1. **IM sim.txt** 每一行为一个场地的模拟，第一列为场地的ID，后面第二列到最后一列为该场地周期（**SiteFile.txt**中的period1）各次随机模拟的结果
2. **IM median with period 0.1.txt** 每一行为一个场地的模拟，第一列为场地的ID，第二列为Sa(T=0.1)的中值地震动强度
3. **IM sim with period 0.1.txt** 每一行为一个场地的模拟，第一列为场地的ID，后面第二列到最后一列为Sa(T=0.1)各次随机模拟的结果

### 示例2 - Python导入

Example 2文件夹包含SpatialIM的Python模块实现。这允许您直接在Python中使用SpatialIM功能。

#### 前提条件

1. **Python 3.12**
2. **相同的处理器架构**：您的Python解释器必须是64位的（win_amd64）。

#### 使用方法

1. 安装:
   ```
   pip install spatialim
   ```
1. 导入:
   ```python
   import spatialim
   ```

2. 创建地震源并设置参数：
   ```python
   # A创建地震源
   eqs = spatialim.EQSource_CB14PCA(lon_0, lat_0)
   
   # 设置参数
   eqs.set_seed(seed)                          # 随机数种子
   eqs.set_W(W)                                # 断裂面宽度(km)
   eqs.set_length(length)                      # 断裂面长度(km)
   eqs.set_RuptureNormal(normal_x, normal_y, normal_z)  # 法线方向
   eqs.set_lambda(lambda_angle)                # 走滑角(度)
   eqs.set_Fhw(Fhw)                            # 上盘效应
   eqs.set_Zhyp(Zhyp)                          # 震源深度(km)
   eqs.set_region(region)                      # 区域(1=加州)
   eqs.set_nPCs(nPCs)                          # 主成分数
   ```

3. 注册场地并运行模拟：
   ```python
   # 注册场地
   eqs.register_site(
       site_id,       # 场地ID
       lon,           # 经度
       lat,           # 纬度
       elevation_km,  # 高程(km)
       period,        # 周期(s)
       vs30,          # Vs30(m/s)
       z25            # Z25(km), 未知用999
   )
   
   # 运行模拟
   magnitudes = [M] * N_sim
   eqs.simulate_intensities(magnitudes, ifmedian=False)
   
   # 保存结果
   eqs.save_im("IM sim.txt")
   eqs.save_xy("XY coord.txt")
   ```

更详细的使用说明，请参考Example 2文件夹中的`README_python_usage.md`和`spatialim_example.py`文件。



## 参考文献

[1] K W Campbell, Y Bozorgnia. NGA-West2 Ground Motion Model for the Average Horizontal Components of PGA, PGV, and 5% Damped Linear Acceleration Response Spectra. Earthquake Spectra, 2014, 30(3): 1087-1115.

[2] N Jayaram, J W Baker. Correlation model for spatially distributed ground-motion intensities. Earthquake Engineering & Structural Dynamics, 2009, 38(15): 1687-1708.

[3] K Goda. Interevent Variability of Spatial Correlation of Peak Ground Motions and Response Spectra. Bulletin of the Seismological Society of America, 2011, 101(5): 2522-2531.

[4] M Markhvida, L Ceferino, J W Baker. Modeling spatially correlated spectral accelerations at multiple periods using principal component analysis and geostatistics. Earthquake Engineering & Structural Dynamics, 2018, 47(5): 1107-1123.
