###############################################################################
# FDTD Plotter for 2D Simulation
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def plot_2d_field(
    filename, 
    Nx, 
    Ny, 
    Lx, 
    Ly, 
    field_name, 
    title="2D FDTD Field"
):
    data = np.loadtxt(filename)
    
    field = data.reshape((Nx, Ny))
    if field_name.startswith('H'):
        field *= 377
    
    x = np.linspace(0, Lx*1e3, Nx)  # Convert to mm 
    y = np.linspace(0, Ly*1e3, Ny)  # Convert to mm
    
    plt.figure(figsize=(10, 8))
    im = plt.imshow(field.T, cmap='RdBu', 
                    aspect='equal',
                    interpolation='nearest',
                    origin='lower',
                    extent=[0, Lx*1e3, 0, Ly*1e3])
    plt.colorbar(im, label=f'{field_name} Magnitude')
    plt.xlabel('X Position (mm)')
    plt.ylabel('Y Position (mm)')
    plt.title(title)
    return plt

def plot_final_frame(directory, Nx, Ny, Lx, Ly, field_name="Array0"):
    """
    Plots only the last snapshot in the directory for the specified field
    
    Parameters
    ----------
    directory : str or Path
        Directory containing the snapshot files
    Nx : int
        Number of grid points in x direction
    Ny : int
        Number of grid points in y direction
    Lx : float
        Physical length in x direction (meters)
    Ly : float
        Physical length in y direction (meters)
    field_name : str
        Name of the field to plot
        
    Returns
    -------
    matplotlib.figure.Figure
        The figure object containing the plot
    """
    # Get all snapshot files for the specified field
    path = Path(directory)
    files = sorted(path.glob(f"{field_name}_*.dat"))
    
    if not files:
        raise FileNotFoundError(f"No snapshot files found for {field_name}")
        
    # Get the last file
    last_file = files[-1]
    frame_num = last_file.stem.split('_')[1]
    
    # Create and save the plot
    plt_obj = plot_2d_field(last_file, Nx, Ny, Lx, Ly, field_name,
                           title=f"{field_name} - Final Frame")
    
    output_dir = path / "plots"
    output_dir.mkdir(exist_ok=True)
    plt_obj.savefig(output_dir / f"{field_name}_final_frame.png")
    
    return plt_obj

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='Plot final frames of FDTD simulation')
    parser.add_argument('--nx', type=int, default=1000,
                        help='Number of grid points in x direction')
    parser.add_argument('--ny', type=int, default=1000, 
                        help='Number of grid points in y direction')
    parser.add_argument('--lx', type=float, default=0.008,
                        help='Physical length in x direction (meters)')
    parser.add_argument('--ly', type=float, default=0.008,
                        help='Physical length in y direction (meters)')
    parser.add_argument('--dir', type=str, default="../data/FDTD2D_1",
                        help='Directory containing snapshot files')
    parser.add_argument('--fields', nargs='+', default=["Ez", "Hy", "Hx"],
                        help='Fields to plot (space separated list)')
    
    args = parser.parse_args()
    
    for field in args.fields:
        plot_final_frame(args.dir, args.nx, args.ny, args.lx, args.ly, field)
