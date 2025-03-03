#include "2DBareBone.h"
#include "snapshot.h"

///============================================================================
/// Simple 2D FDTD simulation with a point source, with a water barrier, 
/// with PML boundary
///============================================================================
int main(int argc, char* argv[]) {
    constexpr std::size_t Nx = 1000;
    constexpr std::size_t Ny = 1000;
    constexpr std::size_t Npmlx = 50;
    constexpr std::size_t Npmly = 50;

    const double Lx = 0.008;
    const double Ly = 0.008;
    const double CourantFactor = 0.99;
    const double tMax = 30e-12;

    const double fSrc = 500e9;
    const double sourcePositionx = 0.004;
    const double sourcePositiony = 0.004;

    PositionType sourcePosition = {sourcePositionx, sourcePositiony};


    std::string dataDir = "../data/2DPointSourcePML"; // Default directory
    if (argc > 1) {
        dataDir = argv[1];
    }

    Snapshot<double> snapshot(dataDir);
    PointSourcePulse<double> source(fSrc, tMax, {sourcePosition});

    BareBone2D<double> simulation(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot
    );

    // Setup PML boundary
    const double m = 10.0;  // PML order
    const double R0 = 1e-8; // Theoretical reflection coefficient
    const double sigmaMax 
        = -(m + 1) * eps0 * c0 * std::log(R0) / (2.0 * std::min(Lx / Nx, Ly / Ny));
    std::cout << "sigmaMax: " << sigmaMax << std::endl;
    simulation.setupPMLBoundary(sigmaMax, m);
    simulation.computeUpdateCoefficients();
    
    std::cout << "Running simulation..." << std::endl;
    simulation.runSimulation();
    std::cout << "Simulation completed." << std::endl;
    return 0;
}