###############################################################################
# FDTD Plotter for 1D Simulation
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def plot_1d_field(filename, field_name, x_points, plt_obj=None, title="1D FDTD Field"):
    data = np.loadtxt(filename)
    
    if field_name.startswith('H'):
        data *= 377
        
    if plt_obj is None:
        plt_obj = plt.figure(figsize=(10, 6))
        plt.xlabel('Distance (mm)')
        plt.ylabel('Field Magnitude')
        plt.title(title)
        plt.grid(True)
        
    x = np.linspace(0, 100, len(data))
    plt.plot(x, data, linewidth=2, label=f'{field_name}')
    plt.legend()
    
    return plt_obj

def plot_final_frame(directory, x_points, field_names=["Array0"]):
    """
    Plots only the last snapshot in the directory for the specified fields,
    overlaying them on the same plot
    
    Parameters
    ----------
    directory : str or Path
        Directory containing the snapshot files
    x_points : int
        Number of grid points in x direction
    field_names : list of str
        Names of the fields to plot
        
    Returns
    -------
    matplotlib.figure.Figure
        The figure object containing the plot
    """
    path = Path(directory)
    plt_obj = None
    
    for field_name in field_names:
        # Get all snapshot files for the specified field
        files = sorted(path.glob(f"{field_name}_*.dat"), 
                      key=lambda x: int(x.stem.split('_')[1]))
        
        if not files:
            raise FileNotFoundError(f"No snapshot files found for {field_name}")
            
        # Get the last file
        last_file = files[-1]
        frame_num = last_file.stem.split('_')[1]
        
        # Create or update the plot
        plt_obj = plot_1d_field(last_file, field_name, x_points, plt_obj,
                               title=f"FDTD Fields - Final Frame ({frame_num})")
    
    output_dir = path / "plots"
    output_dir.mkdir(exist_ok=True)
    plt_obj.savefig(output_dir / "fields_final_frame.png")
    
    return plt_obj


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot final frames of FDTD simulation')
    parser.add_argument('--nx', type=int, default=1000,
                        help='Number of grid points in x direction')
    parser.add_argument('--dir', type=str, default="../data/FDTD1D_1",
                        help='Directory containing snapshot files')
    parser.add_argument('--fields', nargs='+', default=["Ez", "Hy"],
                        help='Fields to plot (space separated list)')
    
    args = parser.parse_args()
    
    plot_final_frame(args.dir, args.nx, args.fields)