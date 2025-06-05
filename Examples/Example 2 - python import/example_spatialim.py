"""
SpatialIM Python示例 - 地震动强度空间分布随机场模拟
===================================================

此示例展示如何使用Python绑定来调用SpatialIM库进行地震动强度空间分布随机场模拟。
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

# ensure the spatialim package has been installed
import spatialim

# 示例：直接使用Python API进行模拟
def example_direct_api():
    print("\n示例：直接使用Python API进行模拟")
    
    # 设置参数
    lon_0, lat_0 = -122.320011, 37.963314  # 震中经纬度
    M = 7.0           # 震级
    N_sim = 100       # 模拟次数
    seed = 42         # 随机数种子
    W = 20.0          # 断裂面宽度
    length = 50.0     # 断裂面长度
    normal_x, normal_y, normal_z = 0.318, 0.214, 0.1395  # 法线方向
    lambda_angle = 0  # rake角度
    Fhw = 1           # hanging wall效应
    Zhyp = 15.0       # 震源深度
    region = 1        # 区域（加州）
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
    return sites

# 主函数
if __name__ == "__main__":
    
    # 运行示例
    sites = example_direct_api()
    
