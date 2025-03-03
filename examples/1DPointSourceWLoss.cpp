#include "1DBareBone.h"
#include "snapshot.h"

///============================================================================
/// Simple 1D FDTD simulation with a point source, but with different material
/// properties in the two regions. In particular, for x < 0.04, the material
/// has a relative permittivity of 4.0, and for x > 0.04, the material has a
/// relative permittivity of 1.0.
///============================================================================
int main(int argc, char* argv[]) {

    /// Simulation parameters
    constexpr std::size_t Nx = 10000;
    constexpr std::size_t Npml = 10;

    const double Lx = 0.10;
    const double CourantFactor = 0.99;
    const double tMax = 100e-12;
    const double fSrc = 30e9;
    const double sourcePosition = 0.05;

    /// Setup simulation data directory
    std::string dataDir = "../data/1DPointSourceWLoss";
    if (argc > 1) {
        dataDir = argv[1];
    }
    Snapshot<double> snapshot(dataDir);

    PointSourcePulse<double> source(fSrc, tMax, {{sourcePosition, 0.0}});
    BareBone1D<double> simulation(Nx, 
                                  Npml,
                                  Lx, 
                                  CourantFactor, 
                                  tMax, 
                                  source, 
                                  snapshot);

    std::vector<double> epsR;
    for (std::size_t i = 0; i < Nx; ++i) {
        if (i < 0.04 * Lx / Nx) {
            epsR.push_back(4.0);
        } else {
            epsR.push_back(1.0);
        }
    }
    simulation.setupMaterialEps(std::move(epsR));
    simulation.computeUpdateCoefficients();
    std::cout << "Running simulation..." << std::endl;
    simulation.runSimulation();
    std::cout << "Simulation completed." << std::endl;
    return 0;
}