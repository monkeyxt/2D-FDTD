#include "2DBareBone.h"
#include "snapshot.h"
#include "pulses.h"
///============================================================================
/// Simple 2D FDTD simulation with a point source to excite TE1 mode in a 
/// planar waveguide, and to study the effect of conductivity on the 
/// propagation of the wave.
///============================================================================
int main(int argc, char* argv[]) {
    constexpr std::size_t Nx = 2000;
    constexpr std::size_t Ny = 2000;
    constexpr std::size_t Npmlx = 50;
    constexpr std::size_t Npmly = 50;

    const double Lx = 0.014;
    const double Ly = 0.014;
    const double CourantFactor = 0.99;
    const double tMax = 8e-11;
    const double fSrc = 175e9;
    const double sourcePositionx = 0.003;
    const double sourcePositiony = 0.008;
    const PositionType sourcePosition = {sourcePositionx, sourcePositiony};


    /// Setup conductivity signature, and fill the regions between x = 2.5 mm 
    /// and x = 14 mm, y = 3.75 mm and y = 4.5 mm, y = 5.5 mm and y = 6.25 mm
    /// with a high conductivity value.
    std::vector<std::vector<double>> sigmaEx(Nx, std::vector<double>(Ny, 0.0));
    
    /// Create the mask for the regions.
    // Convert physical dimensions to grid indices
    const std::size_t x1_lower = static_cast<std::size_t>(0.00650 / Lx * Nx);
    const std::size_t x2_upper = static_cast<std::size_t>(0.01400 / Lx * Nx); 
    const std::size_t y1_lower = static_cast<std::size_t>(0.00375 / Ly * Ny);
    const std::size_t y1_upper = static_cast<std::size_t>(0.00450 / Ly * Ny); 
    const std::size_t y2_lower = static_cast<std::size_t>(0.00550 / Ly * Ny); 
    const std::size_t y2_upper = static_cast<std::size_t>(0.00625 / Ly * Ny); 

    // Run 10 simulations with different conductivity values
    const double sigmaMin = 100.0;
    const double sigmaMax = 1000.0;
    const double sigmaDelta = (sigmaMax - sigmaMin) / 9.0;

    /// Run 10 simulations with different conductivity values
    for (int simIndex = 0; simIndex < 10; simIndex++) {
        const double currentSigma = sigmaMin + simIndex * sigmaDelta;
        
        std::string dataDir = "../data/waveguideAttenuation/sim" + std::to_string(simIndex);
        Snapshot<double> snapshot(dataDir);
        PointSourcePulse<double> source(fSrc, tMax, {sourcePosition});
        BareBone2D<double> simulation(
            Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot
        );

        // Fill the regions with current conductivity value
        std::vector<std::vector<double>> sigmaEx(Nx, std::vector<double>(Ny, 0.0));
        for(std::size_t i = x1_lower; i < x2_upper; i++) {
            for(std::size_t j = y1_lower; j < y1_upper; j++) {
                sigmaEx[i][j] = currentSigma;
            }
            for(std::size_t j = y2_lower; j < y2_upper; j++) {
                sigmaEx[i][j] = currentSigma;
            }
        }
        std::vector<std::vector<double>> sigmaEy(sigmaEx);
        std::vector<std::vector<double>> sigmaMx(sigmaEx);
        std::vector<std::vector<double>> sigmaMy(sigmaEx);
        simulation.setupMaterialSigmaE(std::move(sigmaEx), std::move(sigmaEy));
        simulation.setupMaterialSigmaM(std::move(sigmaMx), std::move(sigmaMy));
        
        /// Setup PML boundary
        const double m = 10.0;  // PML order
        const double R0 = 1e-8; // Theoretical reflection coefficient
        const double pmlSigmaMax 
            = -(m + 1) * eps0 * c0 * std::log(R0) / (2.0 * std::min(Lx / Nx, Ly / Ny));
        std::cout << "PML sigmaMax: " << pmlSigmaMax << std::endl;
        simulation.setupPMLBoundary(pmlSigmaMax, m);
        simulation.computeUpdateCoefficients();
        
        std::cout << "Running simulation " << simIndex << " with sigma = " << currentSigma << "..." << std::endl;
        simulation.runSimulation();
        std::cout << "Simulation " << simIndex << " completed." << std::endl;
    }

    return 0;
}