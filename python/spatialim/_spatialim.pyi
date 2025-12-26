from typing import Optional

def run_simulation(
    eq_source_path: str, 
    site_file_path: str, 
    gmpe_model: Optional[str] = None,
    output_dir: Optional[str] = None
) -> None:
    """
    Run the main simulation workflow. If the number of sites is very large,
    the routine will automatically generate a grid and interpolate results.

    Args:
        eq_source_path (str): Earthquake source file path (.json).
        site_file_path (str): Site file path (.csv).
        gmpe_model (Optional[str], optional): GMPE model name; defaults to "CB14".
        output_dir (Optional[str], optional): Output directory; defaults to "output".
    """
    ...

def run_simulation_grid(
    eq_source_path: str,
    min_lon: float,
    max_lon: float,
    min_lat: float,
    max_lat: float,
    grid_spacing_km: Optional[float] = None,
    gmpe_model: Optional[str] = None,
    output_dir: Optional[str] = None
) -> None:
    """
    Run the simulation in grid mode.
    Generates a grid within the specified longitude/latitude bounds and runs the simulation.

    Args:
        eq_source_path (str): Earthquake source file path (.json).
        min_lon (float): Minimum longitude.
        max_lon (float): Maximum longitude.
        min_lat (float): Minimum latitude.
        max_lat (float): Maximum latitude.
        grid_spacing_km (Optional[float], optional): Grid spacing in kilometers; defaults to 0.5.
        gmpe_model (Optional[str], optional): GMPE model name; defaults to "CB14".
        output_dir (Optional[str], optional): Output directory; defaults to "output".
    """
    ...
