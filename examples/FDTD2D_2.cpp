#include "2DBareBone.h"
#include "snapshot.h"

///============================================================================
/// Simple 2D FDTD simulation with a point source, with a water barrier
///============================================================================
int main(int argc, char* argv[]) {
    constexpr std::size_t Nx = 2000;
    constexpr std::size_t Ny = 2000;
    constexpr std::size_t Npmlx = 50;
    constexpr std::size_t Npmly = 50;

    const double Lx = 0.008;
    const double Ly = 0.008;
    const double CourantFactor = 0.99;
    const double tMax = 10e-12;
    const double fSrc = 500e9;
    const double sourcePositionx = 0.004;
    const double sourcePositiony = 0.004;

    std::string dataDir = "../data/FDTD2D_2"; // Default directory
    if (argc > 1) {
        dataDir = argv[1];
    }

    Snapshot<double> snapshot(dataDir);
    BareBone2D<double> simulation(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, fSrc, 
        sourcePositionx, sourcePositiony, snapshot
    );
    simulation.setupMaterialEps(0.005, 1.0, 1.7);
    simulation.setupMaterialMu(0.005, 1.0, 10.0);
    simulation.computeUpdateCoefficients();
    
    std::cout << "Running simulation..." << std::endl;
    simulation.runSimulation();
    std::cout << "Simulation completed." << std::endl;
    return 0;
}