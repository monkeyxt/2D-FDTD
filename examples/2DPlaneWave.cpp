#include "2DBareBone.h"
#include "snapshot.h"
#include "pulses.h"
///============================================================================
/// Simple 2D FDTD simulation with a plane wave source
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
    const PositionType sourcePosition = {sourcePositionx, sourcePositiony};

    std::string dataDir1 = "../data/2DPlaneWave/1";
    std::string dataDir2 = "../data/2DPlaneWave/2";
    if (argc > 1) {
        dataDir1 = argv[1];
        dataDir2 = argv[2];
    }
    Snapshot<double> snapshot1(dataDir1);
    Snapshot<double> snapshot2(dataDir2);

    /// Initialize a plane wave source
    PlaneWavePulse<double> source(fSrc, tMax, {sourcePosition});

    BareBone2D<double> simulation1(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot1
    );
    simulation1.computeUpdateCoefficients();
    std::cout << "Running simulation 1..." << std::endl;
    simulation1.runSimulation();
    std::cout << "Simulation 1 completed." << std::endl;


    BareBone2D<double> simulation2(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot2
    );

    std::vector<std::vector<double>> epsR(Nx, std::vector<double>(Ny, 1.0));
    for (std::size_t i = 0; i < Nx; i++) {
        for (std::size_t j = 0; j < Ny; j++) {
            if (i * Lx / Nx > 0.005) {
                epsR[i][j] = 3.4;
            }
        }
    }
    simulation2.setupMaterialEps(std::move(epsR));
    simulation2.computeUpdateCoefficients();
    std::cout << "Running simulation 2..." << std::endl;
    simulation2.runSimulation();
    std::cout << "Simulation 2 completed." << std::endl;
    
    return 0;
}