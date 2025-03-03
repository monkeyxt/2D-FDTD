#ifndef EMSOLVER_PULSES_H
#define EMSOLVER_PULSES_H

#include <cmath>
#include <numbers>

#include "constants.h"

///============================================================================
/// Pulse class definition. Defines a pulse with a given frequency, maximum 
/// time, and source position.
///============================================================================
template <FloatingPoint T>
class Pulse {
public:
    Pulse(const T fSrc_, 
          const T tMax_, 
          const std::vector<PositionType<T>> sourcePositions_) 
        : fSrc_(fSrc_), tMax_(tMax_), sourcePositions_(sourcePositions_) {}
    ~Pulse() = default;

    Pulse(const Pulse&) = delete;
    Pulse(Pulse&&) = delete;
    Pulse& operator=(const Pulse&) = delete;
    Pulse& operator=(Pulse&&) = delete;

    /// Compute the pulse at a given time
    /// Base version takes just time parameter
    virtual T computePulse(T t, const PositionTypeIdx& sourceIndex) = 0;

    /// Get the source positions
    std::vector<PositionType<T>> getSourcePositions() const {
        return sourcePositions_;
    }

protected:
    T fSrc_;
    T tMax_;
    std::vector<PositionType<T>> sourcePositions_;
};

///============================================================================
/// Point source pulse class definition
///============================================================================
template <FloatingPoint T>
class PointSourcePulse : public Pulse<T> {
public:
    PointSourcePulse(const T fSrc_, 
                     const T tMax_, 
                     const std::vector<PositionType<T>> sourcePositions_)
        : Pulse<T>(fSrc_, tMax_, sourcePositions_) {}
    ~PointSourcePulse() = default;

    PointSourcePulse(const PointSourcePulse&) = delete;
    PointSourcePulse(PointSourcePulse&&) = delete;
    PointSourcePulse& operator=(const PointSourcePulse&) = delete;
    PointSourcePulse& operator=(PointSourcePulse&&) = delete;

    T computePulse(const T t, const PositionTypeIdx& sourceIndex) override; 
};

template <FloatingPoint T>
T PointSourcePulse<T>::computePulse(const T t, 
                                    const PositionTypeIdx& sourceIndex) {
    return std::sin(2 * std::numbers::pi * this->fSrc_ * t);
}

///============================================================================
/// Gaussian pulse class definition
///============================================================================
template <FloatingPoint T>
class GaussianPointPulse : public PointSourcePulse<T> {
public:
    GaussianPointPulse(const T fSrc_, 
                  const T tMax_, 
                  const std::vector<PositionType<T>> sourcePositions_,
                  const T sigma_ = 1e-12) 
        : Pulse<T>(fSrc_, tMax_, sourcePositions_), sigma_(sigma_) {}
    ~GaussianPointPulse() = default;

    GaussianPointPulse(const GaussianPointPulse&) = delete;
    GaussianPointPulse(GaussianPointPulse&&) = delete;
    GaussianPointPulse& operator=(const GaussianPointPulse&) = delete;
    GaussianPointPulse& operator=(GaussianPointPulse&&) = delete;

    T computePulse(const T t, const PositionTypeIdx& sourceIndex) override;

private:
    T sigma_;
};

template <FloatingPoint T>
T GaussianPointPulse<T>::computePulse(const T t, 
                                      const PositionTypeIdx& sourceIndex) {
    return std::exp(-0.5 * std::pow(t - this->tMax_, 2) / std::pow(this->sigma_, 2));
}


///============================================================================
/// Phase shift pulse class definition
///============================================================================
template <FloatingPoint T>
class PhaseShiftedGaussianPulse : public GaussianPointPulse<T> {
public:
    PhaseShiftedGaussianPulse(const T fSrc_, 
                              const T tMax_, 
                              const std::vector<PositionType<T>> sourcePositions_,
                              const T sigma_,
                              const T incidentAngle_)
        : GaussianPointPulse<T>(fSrc_, tMax_, sourcePositions_, sigma_), 
          incidentAngle_(incidentAngle_) {}
    ~PhaseShiftedGaussianPulse() = default;

    PhaseShiftedGaussianPulse(const PhaseShiftedGaussianPulse&) = delete;
    PhaseShiftedGaussianPulse(PhaseShiftedGaussianPulse&&) = delete;
    PhaseShiftedGaussianPulse& operator=(const PhaseShiftedGaussianPulse&) = delete;
    PhaseShiftedGaussianPulse& operator=(PhaseShiftedGaussianPulse&&) = delete;

    T computePulse(const T t, const PositionTypeIdx& sourceIndex) override;

private:
    T incidentAngle_;
};

template <FloatingPoint T>
T PhaseShiftedGaussianPulse<T>::computePulse(const T t, 
                                             const PositionTypeIdx& sourceIndex) {
    const T timeDelay = (sourceIndex.second * std::sin(this->incidentAngle_) / c0);
    const T tShifted = t - timeDelay;
    return std::exp(-0.5 * std::pow((tShifted - this->tMax_), 2) / std::pow(this->sigma_, 2)) * 
           std::sin(2 * std::numbers::pi * this->fSrc_ * tShifted);
}

///============================================================================
/// Plane wave pulse class definition. For the plane wave pulse, the source 
/// position stores the starting line of the plane wave.
///============================================================================
template <FloatingPoint T>
class PlaneWavePulse : public Pulse<T> {
public:
    PlaneWavePulse(const T fSrc_, 
                   const T tMax_, 
                   const std::vector<PositionType<T>> sourcePositions_) 
        : Pulse<T>(fSrc_, tMax_, sourcePositions_) {}
    ~PlaneWavePulse() = default;

    PlaneWavePulse(const PlaneWavePulse&) = delete;
    PlaneWavePulse(PlaneWavePulse&&) = delete;
    PlaneWavePulse& operator=(const PlaneWavePulse&) = delete;
    PlaneWavePulse& operator=(PlaneWavePulse&&) = delete;

    T computePulse(const T t, const PositionTypeIdx& sourceIndex) override;
};

template <FloatingPoint T>
T PlaneWavePulse<T>::computePulse(const T t, const PositionTypeIdx& sourceIndex) {
    return std::sin(2 * std::numbers::pi * this->fSrc_ * t);
}

///============================================================================
/// Plane wave pulse class definition with phase shift. For the plane wave pulse, 
/// the source position stores the starting line of the plane wave.
///============================================================================
template <FloatingPoint T>
class PlaneWavePulseWithPhaseShift : public PlaneWavePulse<T> {
public:
    PlaneWavePulseWithPhaseShift(const T fSrc_, 
                                 const T tMax_, 
                                 const std::vector<PositionType<T>> sourcePositions_,
                                 const T incidentAngle_)
        : PlaneWavePulse<T>(fSrc_, tMax_, sourcePositions_), 
          incidentAngle_(incidentAngle_) {}
    ~PlaneWavePulseWithPhaseShift() = default;

    PlaneWavePulseWithPhaseShift(const PlaneWavePulseWithPhaseShift&) = delete;
    PlaneWavePulseWithPhaseShift(PlaneWavePulseWithPhaseShift&&) = delete;
    PlaneWavePulseWithPhaseShift& operator=(const PlaneWavePulseWithPhaseShift&) = delete;
    PlaneWavePulseWithPhaseShift& operator=(PlaneWavePulseWithPhaseShift&&) = delete;

    T computePulse(const T t, const PositionTypeIdx& sourceIndex) override;

private:
    T incidentAngle_;
};

template <FloatingPoint T>
T PlaneWavePulseWithPhaseShift<T>::computePulse(const T t, 
                                                const PositionTypeIdx& sourceIndex) {
    const T timeDelay = (sourceIndex.second * std::sin(this->incidentAngle_) / c0);
    const T tShifted = t - timeDelay;
    return std::sin(2 * std::numbers::pi * this->fSrc_ * tShifted);
}

#endif