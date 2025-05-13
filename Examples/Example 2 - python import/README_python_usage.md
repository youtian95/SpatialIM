# SpatialIM Python模块使用指南

## 前提条件

1. **Python 3.12**：您已编译的模块文件名为`spatialim.cp312-win_amd64.pyd`，这表明它是针对Python 3.12编译的。
2. **相同的处理器架构**：您的Python解释器必须是64位的（win_amd64）。

## 导入模块

可以将`.pyd`文件复制到Python脚本所在的目录，然后直接导入：

```python
import spatialim
```

## 使用示例

请参考`spatialim_example.py`文件，了解如何使用SpatialIM模块。

### 基本用法

```python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

# 添加spatialim模块所在的路径，确保可以导入
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

import spatialim

# 设置参数
lon_0, lat_0 = -122.320011, 37.963314  # 震中经纬度
M = 7.0           # 震级
N_sim = 100       # 模拟次数
seed = 42         # 随机数种子
W = 20.0          # 断裂面宽度(km)
length = 50.0     # 断裂面长度(km)
normal_x, normal_y, normal_z = 0.318, 0.214, 0.1395  # 法线方向
lambda_angle = 0  # rake角度(度)
Fhw = 1           # hanging wall效应
Zhyp = 15.0       # 震源深度(km)
region = 1        # 区域(1=加州)
nPCs = 10         # 主成分数
ifmedian = False  # 是否输出中位值

# 定义场地
sites = []
for i in range(10):
    for j in range(10):
        sites.append([
            i*10 + j + 1,                            # ID
            -122.45 + i * 0.01,                      # lon
            37.78 + j * 0.01,                        # lat
            0.0,                                     # elevation_km
            0.1,                                     # T0 (周期)
            500.0,                                   # Vs30
            999.0                                    # Z25 (未知)
        ])

# 创建地震源
eqs = spatialim.EQSource_CB14PCA(lon_0, lat_0)

# 设置随机数种子
eqs.set_seed(seed)

# 设置断裂面参数
eqs.set_W(W)
eqs.set_length(length)
eqs.set_RuptureNormal(normal_x, normal_y, normal_z)
eqs.set_lambda(lambda_angle)
eqs.set_Fhw(Fhw)
eqs.set_Zhyp(Zhyp)
eqs.set_region(region)
eqs.set_nPCs(nPCs)

# 注册场地
for site in sites:
    eqs.register_site(
        site[0],  # ID
        site[1],  # lon
        site[2],  # lat
        site[3],  # elevation_km
        site[4],  # T0
        site[5],  # Vs30
        site[6]   # Z25
    )

# 创建震级列表
magnitudes = [M] * N_sim

# 模拟烈度
print(f"开始模拟，震级={M}，模拟次数={N_sim}...")
eqs.simulate_intensities(magnitudes, ifmedian)

# 获取当前工作目录
work_dir = os.getcwd()

# 保存结果到工作目录
output_file = os.path.join(work_dir, "IM sim.txt")
eqs.save_im(output_file)

# 保存场地坐标到工作目录
coords_file = os.path.join(work_dir, "XY coord.txt")
eqs.save_xy(coords_file)

print(f"模拟完成，结果已保存到 {output_file} 和 {coords_file}")
```

## API参考

请参照`spatialim_bindings.cpp`文件中的Python绑定定义，了解所有可用的类、方法和函数。

### 主要类和方法

- `EQSource_CB14PCA(lon_0, lat_0)` - 创建地震源
- `set_seed(seed)` - 设置随机数种子
- `set_W(W)` - 设置断裂面宽度(km)
- `set_length(length)` - 设置断裂面长度(km)
- `set_RuptureNormal(x, y, z)` - 设置断裂面法线方向
- `set_lambda(lambda_)` - 设置rake角度
- `set_Fhw(Fhw)` - 设置hanging wall效应
- `set_Zhyp(Zhyp)` - 设置震源深度(km)
- `set_region(region)` - 设置区域(1=加州)
- `set_nPCs(nPCs)` - 设置主成分数
- `register_site(ID, lon, lat, elevation_km, T0, Vs30, Z25)` - 注册场地
- `simulate_intensities(magnitudes, ifmedian)` - 模拟烈度
- `save_im(filename)` - 保存烈度结果
- `save_xy(filename)` - 保存场地坐标
