#include "2DBareBone.h"
#include "snapshot.h"
#include "pulses.h"
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

    std::string dataDir = "../data/2DPointSourceWLoss"; // Default directory
    if (argc > 1) {
        dataDir = argv[1];
    }

    Snapshot<double> snapshot(dataDir);

    /// Initialize a point source
    PositionType sourcePosition = {sourcePositionx, sourcePositiony};
    PointSourcePulse<double> source(fSrc, tMax, {sourcePosition});

    BareBone2D<double> simulation(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot
    );

    std::vector<std::vector<double>> eps(Nx, std::vector<double>(Ny, 0.0));
    for (std::size_t i = 0; i < Nx; i++) {
        for (std::size_t j = 0; j < Ny; j++) {
            if (i * Lx / Nx > 0.005) {
                eps[i][j] = 0.04;
            }
        }
    }
    simulation.setupMaterialEps(std::move(eps));
    simulation.computeUpdateCoefficients();
    
    std::cout << "Running simulation..." << std::endl;
    simulation.runSimulation();
    std::cout << "Simulation completed." << std::endl;
    return 0;
}