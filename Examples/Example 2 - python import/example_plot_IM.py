
from pathlib import Path
import spatialim.plotting
import matplotlib.pyplot as plt
import numpy as np

# ensure the spatialim package has been installed
import spatialim

# Set parameters
lon_0, lat_0 = 122.320011, 37.963314  # Epicenter coordinates
M = 7.0           # Magnitude
N_sim = 1       # Number of simulations
seed = 42         # Random seed
W = 20.0          # Fault width
length = 50.0     # Fault length
normal_x, normal_y, normal_z = 1, 0, 1  # Normal direction
lambda_angle = 0  # Rake angle
Fhw = 1           # Whether to consider hanging wall effect
Zhyp = 15.0       # Hypocenter depth
region = 3        # Region
nPCs = 10         # Number of principal components
ifmedian = True  # Whether to output median values only

# Define sites
sites = []
for i in range(100):
    for j in range(100):
        sites.append([
            i*100 + j + 1,                       # ID
            lon_0 + i * 0.01 - 0.5,            # lon
            lat_0 + j * 0.01 - 0.5,            # lat
            0.0,                                # elevation_km
            0.1,                                # T0 (period)
            500.0,                              # Vs30
            999.0                               # Z25 (unknown)
        ])

# Create earthquake source
eqs = spatialim.EQSource_CB14PCA(lon_0, lat_0)

# Set random seed
eqs.set_seed(seed)

# Set fault parameters
eqs.set_W(W)
eqs.set_length(length)
eqs.set_RuptureNormal(normal_x, normal_y, normal_z)
eqs.set_lambda(lambda_angle)
eqs.set_Fhw(Fhw)
eqs.set_Zhyp(Zhyp)
eqs.set_region(region)
eqs.set_nPCs(nPCs)

# Register sites
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

# Create magnitude list
magnitudes = [M] * N_sim

# Simulate intensities
eqs.simulate_intensities(magnitudes, ifmedian)

# Save results
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
# Show plot
plt.show()

