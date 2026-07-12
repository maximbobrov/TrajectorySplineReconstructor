#ifndef SPLINEINTERPOLATOR_H
#define SPLINEINTERPOLATOR_H
#include "globals.h"

// Fits a natural cubic B-spline to one scalar track coordinate (x, y, or z of a
// tracked particle) sampled at successive time steps, then lets the caller
// evaluate the resulting smooth position/velocity/acceleration at any
// (fractional) point along the track. optimizeByGradientDescent() additionally refines
// the spline's control points by gradient descent on a smoothness-penalized
// least-squares objective (objectiveValue()) instead of using the plain interpolating solve.
struct CubicSplineInterpolator
{
    CubicSplineInterpolator(std::vector<double> &vec, std::vector<int> &t);
    ~CubicSplineInterpolator();

    double baseSpline(double x, int derivNum = 0);
    void solveControlPoints();
    void solveControlPointsIterative(size_t itn); // iterative alternative to solveControlPoints(); currently unused
    double positionAt(double n);
    double velocityAt(double n);
    double accelerationAt(double n);
    double evaluateSpline(double x, int derivNum = 0);
    double objectiveValue();
    void optimizeByRandomSearch(size_t itn); // random-search alternative to optimizeByGradientDescent(); currently unused
    int optimizeByGradientDescent(size_t itn);
    void computeObjectiveGradient(double fBefore);

    double m_tolerance;
    double m_optimizationStep;
    size_t m_size;
    std::vector<double> m_vec;
    std::vector<int> m_time;
    double m_dF;
    double m_dt;
    bool m_moreMinNumber; // false if this track is too short (< 4 samples) to fit a spline

    double* m_c;     // fitted spline control points (the actual result), size m_size+3
    double* m_f;     // right-hand side / sample values used while solving, size m_size+3
    double* m_x;     // scratch (RHS during elimination, then solution) for solveControlPoints(), size m_size+3
    double* m_A;     // sub-diagonal of the tridiagonal spline-fit system, size m_size+3
    double* m_B;     // super-diagonal of the tridiagonal spline-fit system, size m_size+3
    double* m_C;     // diagonal of the tridiagonal spline-fit system, size m_size+3
    double* m_cGrad; // per-control-point gradient of objectiveValue(), used by optimizeByGradientDescent()
};

#endif // SPLINEINTERPOLATOR_H