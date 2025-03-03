###############################################################################
# Simple script to compute the attenuation coefficient for a 2D FDTD 
# simulation with a plane wave source.
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse


def wrangle_data(trans_dir, Nx, Ny, field_name="Ez",
                 y_slice_idx=None, x_slice_idx=None):
    """
    Wrangle data from transmission directory
    
    Parameters
    ----------
    trans_dir : str or Path
        Directory containing transmission solution snapshots 
    Nx : int
        Number of grid points in x direction
    Ny : int
        Number of grid points in y direction
    field_name : str, optional
        Name of field to compare (default: "Ez")
    y_slice_idx : int, optional
        Index for y-slice to plot. If None, uses middle of domain
    x_slice_idx : int, optional
        Index for x-slice to plot. If None, uses middle of domain
    """
    trans_path = Path(trans_dir)
    
    # Get final snapshot files for both cases
    trans_files = sorted(trans_path.glob(f"{field_name}_*.dat"),
                        key=lambda x: int(x.stem.split('_')[1]))
    
    if not trans_files:
        raise FileNotFoundError(f"No snapshot files found for {field_name}")
        
    # Load the final frames
    trans_data = np.loadtxt(trans_files[-1])
    
    # Reshape data using provided dimensions
    trans_data = trans_data.reshape((Nx, Ny))
    
    # If no y-slice specified, use middle of domain
    if y_slice_idx is None:
        y_slice_idx = Ny // 2
        
    if x_slice_idx is None:
        x_slice_idx = Nx // 2
        
    # Extract slices
    trans_slice = trans_data[:, y_slice_idx]

    return trans_slice

def compute_attenuation(trans_slice, Nx, Lx=0.008):
    """
    Compute attenuation coefficient from transmission field data
    
    Parameters
    ----------
    trans_slice : numpy.ndarray
        1D slice of transmission field data
    Nx : int
        Number of grid points in x direction
    Lx : float, optional
        Length of domain in meters (default: 0.008 m)
        
    Returns
    -------
    tuple
        Attenuation coefficient (\alpha) in m⁻¹ and fitted parameters (a, \alpha)
    """
    # Create spatial coordinates
    x = np.linspace(0, Lx, Nx)
    
    # Source is in middle of domain
    src_idx = Nx // 2

    # Find peaks in the data
    peaks, _ = find_peaks(trans_slice, height=0.1*np.max(trans_slice))
    
    # Filter for only positive peaks after source position
    valid_peaks = peaks[(trans_slice[peaks] > 0) & (peaks > src_idx)]
    peak_positions = x[valid_peaks]
    peak_values = trans_slice[valid_peaks]

    # Define exponential decay function for fitting
    def exp_decay(x, a, alpha):
        return a * np.exp(-alpha * x)

    # Fit exponential decay to peaks
    popt, _ = curve_fit(exp_decay, peak_positions, peak_values)
    a_fit, alpha_fit = popt

    return alpha_fit, (a_fit, alpha_fit)

def plot_attenuation(trans_slice, Nx, Lx=0.008):
    """
    Plot transmission data and fitted attenuation curve
    
    Parameters
    ----------
    trans_slice : numpy.ndarray
        1D slice of transmission field data
    Nx : int
        Number of grid points in x direction
    Lx : float, optional
        Length of domain in meters (default: 0.008 m)
        
    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the plot
    """
    # Create spatial coordinates
    x = np.linspace(0, Lx, Nx)
    
    # Source is in middle of domain
    src_idx = Nx // 2
    
    # Get attenuation coefficient and fit parameters
    alpha_fit, (a_fit, _) = compute_attenuation(trans_slice, Nx, Lx)
    
    # Define exponential decay function
    def exp_decay(x, a, alpha):
        return a * np.exp(-alpha * x)
    
    # Find peaks for plotting, only after source
    peaks, _ = find_peaks(trans_slice[src_idx:], height=0.1*np.max(trans_slice[src_idx:]))
    peaks = peaks + src_idx  # Adjust indices to account for slicing
    peak_positions = x[peaks]
    peak_values = trans_slice[peaks]

    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(x, trans_slice, 'b-', label='Field data')
    plt.plot(peak_positions, peak_values, 'ro', label='Peaks')
    plt.plot(x[src_idx:], exp_decay(x[src_idx:], a_fit, alpha_fit), 'g--', 
             label=f'Fit: α = {alpha_fit:.2f} m⁻¹')

    plt.xlabel('Position (m)')
    plt.ylabel('Ez')
    plt.title('Field Attenuation Analysis')
    plt.legend()
    plt.grid(True)

    return plt.gcf()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute attenuation coefficient from FDTD solution')
    parser.add_argument('--trans_dir', type=str, required=True,
                        help='Directory containing transmission solution')
    parser.add_argument('--nx', type=int, required=True,
                        help='Number of grid points in x direction')
    parser.add_argument('--ny', type=int, required=True,
                        help='Number of grid points in y direction')
    parser.add_argument('--field', type=str, default="Ez",
                        help='Field to analyze')
    parser.add_argument('--y_slice', type=int, default=None,
                        help='Y-index for slice (default: middle of domain)')
    parser.add_argument('--lx', type=float, default=0.008,
                        help='Length of domain in meters (default: 0.008)')

    args = parser.parse_args()
    
    # Get transmission data using wrangle_data function
    trans_slice = wrangle_data(
        args.trans_dir,
        args.nx,
        args.ny,
        args.field,
        args.y_slice,
    )
    
    # Compute attenuation coefficient
    alpha, (a_fit, _) = compute_attenuation(trans_slice, args.nx, args.lx)
    print(f"\nAttenuation coefficient: {alpha:.2f} m⁻¹")
    path = Path(args.trans_dir)
    output_dir = path / "plots"
    output_dir.mkdir(exist_ok=True)
    plot = plot_attenuation(trans_slice, args.nx, args.lx)
    plot.savefig(output_dir / f"{args.field}_attenuation_analysis.png")
    plt.close()