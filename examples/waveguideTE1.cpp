#include "2DBareBone.h"
#include "snapshot.h"
#include "pulses.h"
///============================================================================
/// Simple 2D FDTD simulation with a point source to excite TE1 mode in a 
/// planar waveguide.
///============================================================================
int main(int argc, char* argv[]) {
    constexpr std::size_t Nx = 2000;
    constexpr std::size_t Ny = 2000;
    constexpr std::size_t Npmlx = 100;
    constexpr std::size_t Npmly = 100;

    const double Lx = 0.014;
    const double Ly = 0.014;
    const double CourantFactor = 0.99;
    const double tMax = 8e-11;
    const double fSrc = 175e9;
    const double sourcePositionx = 0.003;
    const double sourcePositiony = 0.008;
    const PositionType sourcePosition = {sourcePositionx, sourcePositiony};

    std::string dataDir1 = "../data/waveguideTE1/1";
    if (argc > 1) {
        dataDir1 = argv[1];
    }
    Snapshot<double> snapshot1(dataDir1);
    PointSourcePulse<double> source(fSrc, tMax, {sourcePosition});
    BareBone2D<double> simulation1(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot1
    );

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

    // Fill the regions with high conductivity
    for(std::size_t i = x1_lower; i < x2_upper; i++) {
        for(std::size_t j = y1_lower; j < y1_upper; j++) {
            sigmaEx[i][j] = 1e100;
        }
        for(std::size_t j = y2_lower; j < y2_upper; j++) {
            sigmaEx[i][j] = 1e100;
        }
    }
    std::vector<std::vector<double>> sigmaEy(sigmaEx);
    std::vector<std::vector<double>> sigmaMx(sigmaEx);
    std::vector<std::vector<double>> sigmaMy(sigmaEx);
    simulation1.setupMaterialSigmaE(std::move(sigmaEx), std::move(sigmaEy));
    simulation1.setupMaterialSigmaM(std::move(sigmaMx), std::move(sigmaMy));

    /// Setup PML boundary
    const double m = 10.0;  // PML order
    const double R0 = 1e-8; // Theoretical reflection coefficient
    const double sigmaMax 
        = -(m + 1) * eps0 * c0 * std::log(R0) / (2.0 * std::min(Lx / Nx, Ly / Ny));
    std::cout << "sigmaMax: " << sigmaMax << std::endl;
    simulation1.setupPMLBoundary(sigmaMax, m);
    simulation1.computeUpdateCoefficients();
    std::cout << "Running simulation 1..." << std::endl;
    simulation1.runSimulation();
    std::cout << "Simulation 1 completed." << std::endl;

    return 0;
}