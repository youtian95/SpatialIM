
from pathlib import Path
import spatialim.plotting
import matplotlib.pyplot as plt
import numpy as np

# ensure the spatialim package has been installed
import spatialim

# 设置参数
lon_0, lat_0 = 122.320011, 37.963314  # 震中经纬度
M = 7.0           # 震级
N_sim = 1       # 模拟次数
seed = 42         # 随机数种子
W = 20.0          # 断裂面宽度
length = 50.0     # 断裂面长度
normal_x, normal_y, normal_z = 1, 0, 1  # 法线方向
lambda_angle = 0  # rake角度
Fhw = 1           # 是否考虑hanging wall效应
Zhyp = 15.0       # 震源深度
region = 3        # 区域
nPCs = 10         # 主成分数
ifmedian = True  # 是否只输出中位值

# 定义场地
sites = []
for i in range(100):
    for j in range(100):
        sites.append([
            i*100 + j + 1,                       # ID
            lon_0 + i * 0.01 - 0.5,            # lon
            lat_0 + j * 0.01 - 0.5,            # lat
            0.0,                                # elevation_km
            0.1,                                # T0 (周期)
            500.0,                              # Vs30
            999.0                               # Z25 (未知)
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
eqs.simulate_intensities(magnitudes, ifmedian)

# 保存结果
work_dir = Path(__file__).parent.resolve()
output_file = work_dir / "IM sim.txt"
eqs.save_im(str(output_file))
coords_file = work_dir / "XY coord.txt"
eqs.save_xy(str(coords_file))


plotter = spatialim.plotting.SpatialIMPlotter()

fig = plotter.plot_intensity_contour(
    str(output_file),
    str(coords_file),
    fault_width=W,
    fault_length=length,
    rupture_normal=np.array([normal_x, normal_y, normal_z]),
    show_fault=True
)
# 显示图形
plt.show()

