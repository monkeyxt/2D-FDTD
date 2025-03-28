#ifndef EMSOLVER_1DBAREBONE_H
#define EMSOLVER_1DBAREBONE_H

#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <functional>
#include <numbers>
#include <cassert>
#include <omp.h>

#include "constants.h"
#include "snapshot.h"
#include "pulses.h"

///============================================================================
/// 1D Barebone Solver for the 1D Wave Equation
///============================================================================
template <FloatingPoint T>
class BareBone1D {
public:
    /// Constructor
    explicit BareBone1D(std::size_t Nx_,
                        std::size_t Npml_,
                        T Lx_, 
                        T courantFactor_, 
                        T tMax_,
                        Pulse<T>& source_,
                        Snapshot<T>& snapshot_);

    ~BareBone1D() = default;

    BareBone1D() = delete;
    BareBone1D(const BareBone1D&) = delete;
    BareBone1D(BareBone1D&&) = delete;
    BareBone1D& operator=(const BareBone1D&) = delete;
    BareBone1D& operator=(BareBone1D&&) = delete;

    /// Setup electric permittivity
    void setupMaterialEps(const std::vector<T>&& epsR_);

    /// Setup magnetic permeability
    void setupMaterialMu(const std::vector<T>&& muR_);

    /// Setup PML boundary
    void setupPMLBoundary(const T sigmaMax, const T factor);

    /// Compute update coefficients
    void computeUpdateCoefficients();

    /// Run simulation
    void runSimulation();

    /// Write fields to snapshot
    void writeFields();

private:
    /// Simulation parameters
    std::size_t Nx;    /// number of grid points
    std::size_t Npml;  /// number of PML points
    T Lx;              /// total length in meters
    T courantFactor;   /// Courant factor
    T dt;              /// time step in seconds
    T tMax;            /// maximum time in seconds
    T dx;              /// grid spacing in meters
    std::size_t Nt;    /// number of time steps

    /// Source parameters
    Pulse<T>& source;  /// source pulse
    std::size_t srcIndex;   /// source index

    /// Field arrays
    std::vector<T> Ez;
    std::vector<T> Hy;

    /// Material parameters
    std::vector<T> epsR;
    std::vector<T> muR;
    std::vector<T> sigmaE;
    std::vector<T> sigmaM;

    /// Update coefficients
    std::vector<T> ceze;
    std::vector<T> cezh;
    std::vector<T> chye;
    std::vector<T> chyh;

    /// Snapshot handler
    Snapshot<T>& snapshot;

    /// Initialization
    void initialize();

    /// Update source
    void updateSource(const std::size_t t);
};

template <FloatingPoint T>
BareBone1D<T>::BareBone1D(std::size_t Nx_,
                          std::size_t Npml_,
                          T Lx_, 
                          T courantFactor_, 
                          T tMax_, 
                          Pulse<T>& source_,
                          Snapshot<T>& snapshot_)
    : Nx(Nx_), Npml(Npml_), Lx(Lx_), courantFactor(courantFactor_), 
      tMax(tMax_), source(source_), snapshot(snapshot_) {
    if (Nx <= 2) {
        throw std::runtime_error("Nx must be greater than 2");
    }
    if (Npml <= 0) {
        throw std::runtime_error("Npml must be greater than 0");
    }
    if (Nx <= Npml) {
        throw std::runtime_error("Nx must be greater than Npml");
    }
    initialize();
}

template <FloatingPoint T>
void BareBone1D<T>::initialize() {
    dx = Lx / Nx;
    dt = courantFactor * dx / c0;
    Nt = static_cast<std::size_t>(tMax / dt);

    /// Zero initial fields
    Ez.resize(Nx);
    Hy.resize(Nx);

    /// Setup material parameters
    epsR.resize(Nx, 1.0);
    muR.resize(Nx, 1.0);
    
    /// Setup conductivity parameters
    sigmaE.resize(Nx);
    sigmaM.resize(Nx);

    /// Setup update coefficients
    ceze.resize(Nx);
    cezh.resize(Nx);
    chye.resize(Nx);
    chyh.resize(Nx);

    /// Only point source is supported for now
    if (auto* pointSource = dynamic_cast<PointSourcePulse<T>*>(&source)) {
        srcIndex = static_cast<std::size_t>(
            pointSource->getSourcePositions()[0].first / dx);
    }
    else {
        throw std::runtime_error("Only point source is supported for now");
    }

    /// Calculate the source index
    srcIndex = static_cast<std::size_t>(sourcePosition / dx);
    if (srcIndex >= Nx) {
        throw std::runtime_error("Source position is out of bounds");
    }

    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Lx = " << Lx << " m" << std::endl;
    std::cout << "dx = " << dx << " m" << std::endl;
    std::cout << "dt = " << dt << " s" << std::endl;
    std::cout << "steps= " << Nt << std::endl;
    std::cout << "src= " << sourcePosition << " m" << std::endl;
    std::cout << "srcIndex= " << srcIndex << std::endl;
    std::cout << "--------------------------------" << std::endl;
}

template <FloatingPoint T>
void BareBone1D<T>::setupMaterialEps(const std::vector<T>&& epsR_) {
    epsR = epsR_;
}

template <FloatingPoint T>
void BareBone1D<T>::setupMaterialMu(const std::vector<T>&& muR_) {
    muR = muR_;
}

template <FloatingPoint T>
void BareBone1D<T>::setupPMLBoundary(const T sigmaMax, const T factor) {
    // Left PML region
    for(std::size_t i = 0; i < Npml; i++) {
        sigmaE[i] = sigmaMax * std::pow((Npml - static_cast<T>(i)) / Npml, factor);
        sigmaM[i] = mu0/eps0 * sigmaE[i];
    }

    // Right PML region
    for(std::size_t i = Nx - Npml; i < Nx; i++) {
        sigmaE[i] = sigmaMax * std::pow(static_cast<T>(i - (Nx - Npml)) / Npml, factor);
        sigmaM[i] = mu0/eps0 * sigmaE[i];
    }
}

template <FloatingPoint T>
void BareBone1D<T>::computeUpdateCoefficients() {
    /// Conduction from Schneider 3.52
    for(std::size_t i = 0; i < Nx; i++) {
        T denomE = 1.0 + (sigmaE[i] * dt) / (2.0 * epsR[i] * eps0);
        T numerE = 1.0 - (sigmaE[i] * dt) / (2.0 * epsR[i] * eps0);
        ceze[i] = numerE / denomE;
        cezh[i] = (dt / (epsR[i] * eps0 * dx)) / denomE;
    }
    /// Conduction from Schneider 3.56
    for(std::size_t i = 0; i < Nx; i++) {
        T denomH = 1.0 + (sigmaM[i] * dt) / (2.0 * muR[i] * mu0);
        T numerH = 1.0 - (sigmaM[i] * dt) / (2.0 * muR[i] * mu0);
        chye[i] = numerH / denomH;
        chyh[i] = (dt / (muR[i] * mu0 * dx)) / denomH;
    }
}

template <FloatingPoint T>
void BareBone1D<T>::runSimulation() {
    /// Run simulation for Nt time steps
    for (std::size_t t = 0; t < Nt; t++) {
        #pragma omp simd
        for (std::size_t i = 0; i < Nx - 1; i++) {
            Hy[i] = chye[i] * Hy[i] + chyh[i] * (Ez[i + 1] - Ez[i]);
        }

        #pragma omp simd
        for (std::size_t i = 1; i < Nx; i++) {
            Ez[i] = ceze[i] * Ez[i] + cezh[i] * (Hy[i] - Hy[i - 1]);
        }
        updateSource(t);
    }
    writeFields();
}

template <FloatingPoint T>
void BareBone1D<T>::updateSource(const std::size_t t) {
    for(const auto& src : srcIndex) {
        auto pulseInc = source.computePulse(t * dt, src);
        Ez[src.first] += pulseInc;
    }
}

template <FloatingPoint T>
void BareBone1D<T>::writeFields() {
    snapshot.write({"Ez", "Hy"}, Ez, Hy);
}

#endif // EMSOLVER_1DBAREBONE_H