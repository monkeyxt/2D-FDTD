#include "constants.h"
#include "2DBareBone.h"
#include "snapshot.h"
#include "pulses.h"

///============================================================================
/// Simple 2D FDTD simulation with a point source
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
    
    PositionType sourcePosition = {sourcePositionx, sourcePositiony};
    

    std::string dataDir = "../data/2DPointSource";
    if (argc > 1) {
        dataDir = argv[1];
    }
    Snapshot<double> snapshot(dataDir);

    // Construct a plane wave pulse
    PointSourcePulse<double> source(fSrc, tMax, {sourcePosition});

    BareBone2D<double> simulation(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot
    );
    simulation.computeUpdateCoefficients();
    std::cout << "Running simulation..." << std::endl;
    simulation.runSimulation();
    std::cout << "Simulation completed." << std::endl;
    return 0;
}