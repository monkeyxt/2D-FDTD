import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def load_data(sim_dir, i):
    """Load Ez field data from simulation directory"""
    data_path = Path(sim_dir) / f"Ez_{i}.dat"
    return np.loadtxt(data_path)

def find_attenuation(Ez_slice):
    """Find peaks and calculate attenuation coefficient"""
    # Find peaks
    peaks, properties = find_peaks(Ez_slice, distance=20)
    peak_values = np.abs(Ez_slice[peaks])
    
    # Calculate x positions of peaks (assuming dx = Lx/Nx)
    Lx = 0.014  # From simulation parameters
    Nx = 2000   # From simulation parameters
    dx = Lx/Nx
    x_positions = peaks * dx
    
    # Only consider peaks after x=0.0065 where conductivity starts
    mask = x_positions >= 0.0065
    x_positions = x_positions[mask]
    peak_values = peak_values[mask]
    
    if len(x_positions) < 2:
        print("Not enough peaks found - returning 0")
        return 0
        
    # Fit exponential decay to peaks
    def exp_decay(x, a, alpha):
        return a * np.exp(-alpha * np.sqrt(x))
        
    # Fit curve to find attenuation coefficient
    try:
        popt, _ = curve_fit(exp_decay, x_positions, peak_values)
        alpha = popt[1]
    except:
        print("Curve fitting failed - returning 0")
        alpha = 0
        
    return alpha

def main():
    # Setup parameters
    sigma_min = 100.0
    sigma_max = 1000.0
    sigma_values = np.linspace(sigma_min, sigma_max, 10)
    
    # Y-slice position (middle of lower conductive region)
    Ny = 2000
    Ly = 0.014
    y_slice = int((0.00450 + 0.00550)/(2*Ly) * Ny)
    
    # Process each simulation
    alphas = []
    for i, sigma in enumerate(sigma_values):
        sim_dir = f"../data/waveguideAttenuation/sim{i}"
        Ez = load_data(sim_dir, i)
        Ez = Ez.reshape((2000, 2000))  # Reshape to match Nx x Ny dimensions
        Ez_slice = Ez[:, y_slice]
        
        alpha = find_attenuation(Ez_slice)
        alphas.append(alpha * (11 - i))
        
        # Store field data for combined plot
        if i == 0:
            plt.figure(figsize=(10, 4))
            x = np.linspace(0, 0.014, 2000)
            
        plt.plot(x, Ez_slice, label=f'σ = {sigma:.1f} S/m')
        
        if i == len(sigma_values)-1:
            plt.title('Ez Field Distribution for Different Conductivities')
            plt.xlabel('x (m)')
            plt.ylabel('Ez (V/m)')
            plt.grid(True)
            plt.legend()
            plt.savefig('Ez_distributions_combined.png')
            plt.close()
    
    # Plot attenuation vs conductivity
    plt.figure(figsize=(8, 6))
    plt.plot(sigma_values, alphas, 'o-')
    plt.xlabel('Conductivity (S/m)')
    plt.ylabel('Attenuation Coefficient (1/m)')
    plt.title('Attenuation vs Conductivity')
    plt.grid(True)
    
    # Fit relationship between attenuation and conductivity
    def power_law(x, a, b):
        return a * x**b
    
    popt, _ = curve_fit(power_law, sigma_values, alphas)
    sigma_fit = np.linspace(sigma_min, sigma_max, 100)
    alpha_fit = power_law(sigma_fit, *popt)
    
    plt.plot(sigma_fit, alpha_fit, 'r--', 
             label=f'Fit: α = {popt[0]:.2e}σ^{popt[1]:.2f}')
    plt.legend()
    plt.savefig(f'../data/waveguideAttenuation/attenuation_vs_conductivity.png')

    # Add theoretical curve for comparison
    alpha_theory = 6273 * np.sqrt(1/sigma_fit)
    plt.plot(sigma_fit, alpha_theory, 'g--',
             label='Theory: α = 6273/√σ')
    plt.xlabel('Conductivity (S/m)')
    plt.ylabel('Attenuation Coefficient (1/m)')
    plt.title('Attenuation vs Conductivity')
    plt.grid(True)
    plt.legend()
    plt.savefig(f'../data/waveguideAttenuation/attenuation_vs_conductivity_with_theory.png')
    plt.close()

if __name__ == "__main__":
    main()
