# Programs for Simulating Ground Motion Intensity Measure Random Fields

[中文文档](./README_CN.md)

## Examples
Below are the spectral acceleration (Sa) distributions at T=0.2, 0.5, and 1.0s for a region under an earthquake:
![Sa0.2](./Figures/Sa0.2.png)
![Sa0.5](./Figures/Sa0.5.png)
![Sa1.0](./Figures/Sa1.0.png)

## Usage

### Examples 1 - executable program

The Examples 1 folder contains sample files. Place two files named `EQSource.txt` and `SiteFile.txt`, containing earthquake source information and site information respectively, and then run the `IMSim.exe` program to perform the simulation.

#### Input

1. **EQSource.txt** Each line contains the following parameters (separated by spaces):
    - ifmedian - 0/1, whether to output median values
    - M - Magnitude
    - N_sim - Number of simulations
    - seed - int, random number seed
    - lon_0, lat_0 - Epicenter longitude and latitude, degrees
    - W - Width of the fault rupture plane rectangle, km, input 999 if unknown
    - length - Length of the fault rupture plane rectangle, km
    - RuptureNormal_x, RuptureNormal_y, RuptureNormal_z - Normal direction of fault plane (east is x, north is y, up is z)
    - lambda - Rake angle (°) - Average slip angle measured on the hanging wall within the rupture plane, 0 degrees aligned with strike direction, positive counterclockwise
    - Fhw - Whether to consider hanging wall effects, 0/1
    - Zhyp - Hypocenter depth measured from sea level, km, input 999 if unknown
    - region - Study region
     = 0 Global (including Taiwan)
     = 1 California
     = 3 China or Turkey
     = 4 Italy
    - nPCs - Number of principal components to consider in the IM correlation PCA method, recommended 5 or more
1. **SiteFile.txt** Each line represents data for one site, with the following parameters (separated by spaces):
    - ID - Site ID
    - lon - Longitude
    - lat - Latitude
    - elevation_km - Elevation, km
    - period1 - Fundamental period
    - Vs30_mpers - Shear wave velocity
    - Z25_km - Depth to the 2.5km/s shear wave velocity horizon, km (if in California or Japan and Z25_km is unknown, input 999)

#### Output

1. **IM sim.txt** 
   Each line represents simulation results for one site. The first column is the site ID, and subsequent columns (from the second to the last) are the random simulation results for the site's period (period1 in `SiteFile.txt`)
2. **IM median with period 0.1.txt** 
   Each line represents simulation results for one site. The first column is the site ID, and the second column is the median seismic intensity for $Sa(T=0.1)$
3. **IM sim with period 0.1.txt** 
   Each line represents simulation results for one site. The first column is the site ID, and subsequent columns (from the second to the last) are the random simulation results for $Sa(T=0.1)$

### Example 2 - Python Import

The Example 2 folder contains a Python module implementation of SpatialIM. This allows you to use the SpatialIM functionality directly within Python.

#### Requirements

1. **Python 3.12**
2. **Matching Processor Architecture**: Your Python interpreter must be 64-bit (win_amd64).

#### Usage
1. Install it:
   ```
   pip install spatialim
   ```
1. Import it:
   ```python
   import spatialim
   ```

1. Create an earthquake source and set parameters:
   ```python
   # Create earthquake source
   eqs = spatialim.EQSource_CB14PCA(lon_0, lat_0)
   
   # Set parameters
   eqs.set_seed(seed)                          # Random number seed
   eqs.set_W(W)                                # Fault rupture width (km)
   eqs.set_length(length)                      # Fault rupture length (km)
   eqs.set_RuptureNormal(normal_x, normal_y, normal_z)  # Normal direction
   eqs.set_lambda(lambda_angle)                # Rake angle (degrees)
   eqs.set_Fhw(Fhw)                            # Hanging wall effect
   eqs.set_Zhyp(Zhyp)                          # Hypocenter depth (km)
   eqs.set_region(region)                      # Region (1=California)
   eqs.set_nPCs(nPCs)                          # Number of principal components
   ```

5. Register sites and run simulations:
   ```python
   # Register a site
   eqs.register_site(
       site_id,       # ID
       lon,           # Longitude
       lat,           # Latitude
       elevation_km,  # Elevation (km)
       period,        # Period (s)
       vs30,          # Vs30 (m/s)
       z25            # Z25 (km), use 999 if unknown
   )
   
   # Run simulations
   magnitudes = [M] * N_sim
   eqs.simulate_intensities(magnitudes, ifmedian=False)
   
   # Save results
   eqs.save_im("IM sim.txt")
   eqs.save_xy("XY coord.txt")
   ```

For more detailed instructions, refer to the `README_python_usage.md` and `spatialim_example.py` files in the Example 2 folder.

## References

[1] K W Campbell, Y Bozorgnia. NGA-West2 Ground Motion Model for the Average Horizontal Components of PGA, PGV, and 5% Damped Linear Acceleration Response Spectra. Earthquake Spectra, 2014, 30(3): 1087-1115.

[2] N Jayaram, J W Baker. Correlation model for spatially distributed ground-motion intensities. Earthquake Engineering & Structural Dynamics, 2009, 38(15): 1687-1708.

[3] K Goda. Interevent Variability of Spatial Correlation of Peak Ground Motions and Response Spectra. Bulletin of the Seismological Society of America, 2011, 101(5): 2522-2531.

[4] M Markhvida, L Ceferino, J W Baker. Modeling spatially correlated spectral accelerations at multiple periods using principal component analysis and geostatistics. Earthquake Engineering & Structural Dynamics, 2018, 47(5): 1107-1123.
