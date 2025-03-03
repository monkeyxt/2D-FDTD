###############################################################################
# Simple script to compute the transmission coefficient for a 2D FDTD 
# simulation with a plane wave source.
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def wrangle_data(ref_dir, trans_dir, Nx, Ny, field_name="Ez",
                 y_slice_idx=None, x_slice_idx=None):
    """
    Wrangle data from reference and transmission directories
    
    Parameters
    ----------
    ref_dir : str or Path
        Directory containing reference solution snapshots
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
    ref_path = Path(ref_dir)
    trans_path = Path(trans_dir)
    
    # Get final snapshot files for both cases
    ref_files = sorted(ref_path.glob(f"{field_name}_*.dat"), 
                      key=lambda x: int(x.stem.split('_')[1]))
    trans_files = sorted(trans_path.glob(f"{field_name}_*.dat"),
                        key=lambda x: int(x.stem.split('_')[1]))
    
    if not ref_files or not trans_files:
        raise FileNotFoundError(f"No snapshot files found for {field_name}")
        
    # Load the final frames
    ref_data = np.loadtxt(ref_files[-1])
    trans_data = np.loadtxt(trans_files[-1])
    
    # Reshape data using provided dimensions
    ref_data = ref_data.reshape((Nx, Ny))
    trans_data = trans_data.reshape((Nx, Ny))
    
    # If no y-slice specified, use middle of domain
    if y_slice_idx is None:
        y_slice_idx = Ny // 2
        
    if x_slice_idx is None:
        x_slice_idx = Nx // 2
        
    # Extract slices
    ref_slice = ref_data[:, y_slice_idx]
    trans_slice = trans_data[:, y_slice_idx]

    return ref_slice, trans_slice
    
def plot_field_comparison(ref_slice, trans_slice, Nx, field_name="Ez"):
    # Create x coordinates in mm (assuming Lx = 8mm from example files)
    Lx = 8.0  # Length in mm
    x = np.linspace(0, Lx, Nx)
    
    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(x, ref_slice, 'b-', label='Reference')
    plt.plot(x, trans_slice, 'r--', label='Transmission')
    plt.xlabel('x (mm)')
    plt.ylabel(f'{field_name} Field Amplitude')
    plt.title(f'{field_name} Field Comparison Along y-slice')
    plt.legend()
    plt.grid(True)
    
    return plt.gcf()

def compute_transmission_coefficient(ref_slice, trans_slice, source_idx=None):
    """
    Compute the transmission coefficient from reference and transmission field data
    
    Parameters
    ----------
    ref_slice : numpy.ndarray
        1D slice of reference field data
    trans_slice : numpy.ndarray
        1D slice of transmission field data
    source_idx : int, optional
        Index of source location. If None, uses 1/3 of domain length
        
    Returns
    -------
    float
        Transmission coefficient (|T|)
    """
    if source_idx is None:
        source_idx = len(ref_slice) // 3
        
    incident_amp = np.max(np.abs(ref_slice[source_idx:source_idx+100]))
    transmitted_amp = np.max(np.abs(trans_slice[source_idx:source_idx+100]))
    return np.abs(transmitted_amp / incident_amp)

def compute_reflection_coefficient(ref_slice, trans_slice, source_idx=None):
    """
    Compute the reflection coefficient from reference and transmission field data
    
    Parameters
    ----------
    ref_slice : numpy.ndarray
        1D slice of reference field data
    trans_slice : numpy.ndarray
        1D slice of transmission field data
    source_idx : int, optional
        Index of source location. If None, uses 1/3 of domain length
        
    Returns
    -------
    float
        Reflection coefficient (|Γ|)
    """
    if source_idx is None:
        source_idx = len(ref_slice) // 3
        
    incident_amp = np.max(np.abs(ref_slice[0:source_idx-100]))
    reflected_amp = np.max(np.abs(trans_slice[0:source_idx-100] - ref_slice[0:source_idx-100]))
    return np.abs(reflected_amp / incident_amp)

def compute_standing_wave_ratio(ref_slice, trans_slice, source_idx=None):
    """
    Compute the standing wave ratio (SWR) from the field data
    
    Parameters
    ----------
    ref_slice : numpy.ndarray
        1D slice of reference field data
    trans_slice : numpy.ndarray
        1D slice of transmission field data
    source_idx : int, optional
        Index of source location. If None, uses 1/3 of domain length
        
    Returns
    -------
    float
        Standing Wave Ratio (SWR)
    """
    reflection_coef = compute_reflection_coefficient(ref_slice, trans_slice, source_idx)
    return (1 + reflection_coef) / (1 - reflection_coef)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare transmission between two FDTD solutions')
    parser.add_argument('--ref_dir', type=str, required=True,
                        help='Directory containing reference solution')
    parser.add_argument('--trans_dir', type=str, required=True,
                        help='Directory containing transmission solution')
    parser.add_argument('--nx', type=int, required=True,
                        help='Number of grid points in x direction', default=2000)
    parser.add_argument('--ny', type=int, required=True,
                        help='Number of grid points in y direction', default=2000)
    parser.add_argument('--field', type=str, default="Ez",
                        help='Field to compare')
    parser.add_argument('--y_slice', type=int, default=None,
                        help='Y-index for slice (default: middle of domain)')
    parser.add_argument('--x_slice', type=int, default=None,
                        help='X-index for slice (default: middle of domain)')

    args = parser.parse_args()
    
    ref_slice, trans_slice = wrangle_data(
        args.ref_dir, 
        args.trans_dir, 
        args.nx, 
        args.ny, 
        args.field, 
        args.y_slice,
        args.x_slice
    )
    
    plot = plot_field_comparison(ref_slice, trans_slice, args.nx, args.field)
    path = Path(args.ref_dir)
    output_dir = path / "plots"
    output_dir.mkdir(exist_ok=True)
    plot.savefig(output_dir / f"{args.field}_transmission_comparison.png")
    
    src_idx = len(ref_slice) // 2
    trans_coef = compute_transmission_coefficient(ref_slice, trans_slice, 
                                                  src_idx)
    refl_coef = compute_reflection_coefficient(ref_slice, trans_slice, 
                                               src_idx)
    swr = compute_standing_wave_ratio(ref_slice, trans_slice, 
                                      src_idx)
    
    print(f"Transmission Coefficient: {trans_coef:.4f}")
    print(f"Reflection Coefficient: {refl_coef:.4f}")
    print(f"Standing Wave Ratio: {swr:.4f}")