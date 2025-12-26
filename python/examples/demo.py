import spatialim
from spatialim.plotting import plot_intensity_contour
import os
import json
import pandas as pd
import numpy as np

def main():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    eq_source_path = os.path.join(current_dir, "fixtures", "EQSource.json")
    output_dir = os.path.join(current_dir, "output")
    
    spatialim.run_simulation_grid(
        eq_source_path=eq_source_path,
        min_lon=99.5, max_lon=100.5,
        min_lat=29.5, max_lat=30.5,
        grid_spacing_km=1.0,
        output_dir=output_dir
    )
    
    # Read the actual source parameters used by the Rust core
    used_eq_source_path = os.path.join(output_dir, "used_eq_source.json")
    with open(used_eq_source_path, 'r') as f:
        eq_data = json.load(f)
            

    # Plot the T = 1 s results
    im_file = os.path.join(output_dir, "IM_T1.csv")
    if os.path.exists(im_file):
        df = pd.read_csv(im_file)
        
        # 1) Plot the first simulation (Sim_1)
        # CSV format is assumed: Site_ID, Sim_1, Sim_2, ...
        # Keep only Site_ID and Sim_1 for the first plot
        temp_im_file_sim1 = os.path.join(output_dir, "temp_im_sim1.txt")
        df.iloc[:, [0, 1]].to_csv(temp_im_file_sim1, sep=' ', index=False, header=False)
        
        plot_intensity_contour(
            im_file=temp_im_file_sim1,
            coords_file=os.path.join(output_dir, "grid_coordinates.csv"),
            title="PSA (T=1 s) Simulation 1",
            save_path=os.path.join(output_dir, "plot_sim1.png"),
            fault_strike=eq_data.get('strike'),
            fault_dip=eq_data.get('dip'),
            fault_width=eq_data.get('w'),
            fault_length=eq_data.get('length'),
            show_fault=True
        )

        # 2) Plot the median over all simulations (log-mean, then exp)
        if df.shape[1] > 1:
            # Compute log-mean and exponentiate back
            sim_cols = df.iloc[:, 1:]
            median_values = np.exp(np.log(sim_cols).mean(axis=1))
            
            # Build DataFrame for plotting: Site_ID, Median
            df_median = pd.DataFrame({'Site_ID': df.iloc[:, 0], 'Median': median_values})
            
            temp_im_file_median = os.path.join(output_dir, "temp_im_median.txt")
            df_median.to_csv(temp_im_file_median, sep=' ', index=False, header=False)
            
            plot_intensity_contour(
                im_file=temp_im_file_median,
                coords_file=os.path.join(output_dir, "grid_coordinates.csv"),
                title="PSA (T=1 s) Median",
                save_path=os.path.join(output_dir, "plot_median.png"),
                fault_strike=eq_data.get('strike'),
                fault_dip=eq_data.get('dip'),
                fault_width=eq_data.get('w'),
                fault_length=eq_data.get('length'),
                show_fault=True
            )


if __name__ == "__main__":
    main()
