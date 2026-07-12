#include "splineinterpolator.h"

CubicSplineInterpolator::CubicSplineInterpolator(std::vector<double> &vec, std::vector<int> &t)
{
    m_size = vec.size();
    m_moreMinNumber = vec.size() >= 4;
    m_vec = vec;
    m_time.resize(t.size());
    for (size_t i = 0; i < t.size(); i++)
    {
        m_time[i] = t[i];
    }
    m_c     = new double[m_size + 3];
    m_f     = new double[m_size + 3];
    m_x     = new double[m_size + 3];
    m_A     = new double[m_size + 3];
    m_B     = new double[m_size + 3];
    m_C     = new double[m_size + 3];
    m_cGrad = new double[m_size + 3];

    for (size_t i = 1; i < m_size + 1; i++)
    {
        m_f[i] = vec[i - 1];
        m_x[i] = vec[i - 1];
    }
    m_dF = 1.0;
    m_tolerance = 1e-6;
    m_dt = 0.0006;
    m_optimizationStep = m_tolerance;
    solveControlPoints();
    // scratch arrays only needed while solving; m_c (the fitted result) and
    // m_cGrad (used later by optimizeByGradientDescent()) are kept.
    delete[] m_f;
    delete[] m_x;
    delete[] m_A;
    delete[] m_B;
    delete[] m_C;
}

CubicSplineInterpolator::~CubicSplineInterpolator()
{
    delete[] m_c;
    delete[] m_cGrad;
}

// Cubic B-spline basis function and its derivatives (up to 2nd order),
// evaluated on its compact support x in [-2, 2].
double CubicSplineInterpolator::baseSpline(double x, int derivNum )
{
    if (derivNum == 0)
    {
        if (x <= -2.0)
        {
            return 0;
        }
        else if ((x > -2.0) && (x <= -1.0))
        {
            return (1.0 / 6.0) * (2.0 + x) * (2.0 + x) * (2.0 + x);
        }
        else if ((x > -1.0) && (x <= 0.0))
        {
            return (4.0 / 6.0 - (x) * (x) -0.5 * (x) * (x) * (x));
        }
        else if ((x > 0.0) && (x <= 1.0))
        {
            return (4.0 / 6.0 - x * x + 0.5 * x * x * x);
        }
        else if ((x > 1.0) && (x <= 2.0))
        {
            return (1.0 / 6.0) * (2.0 - x) * (2.0 - x) * (2.0 - x);
        }
        else if (x > 2.0)
        {
            return 0;
        }
    }
    if (derivNum == 1)
    {
        if (x <= -2.0)
        {
            return 0;
        }
        else if ((x > -2.0) && (x <= -1.0))
        {
            return (3.0 / 6.0) * (2.0 + x) * (2.0 + x);
        }
        else if ((x > -1.0) && (x <= 0))
        {
            return (-2.0 * x - 1.5 * x * x);
        }
        else if ((x > 0) && (x <= 1.0))
        {
            return (-2.0 * x + 1.5 * x * x);
        }
        else if ((x > 1.0) && (x <= 2.0))
        {
            return (-3.0 / 6.0) * (2.0 - x) * (2.0 - x);
        }
        else if (x > 2.0)
        {
            return 0;
        }
    }
    if (derivNum == 2)
    {
        if (x <= -2.0)
        {
            return 0;
        }
        else if ((x > -2.0) && (x <= -1.0))
        {
            return (2.0 + x);
        }
        else if ((x > -1.0) && (x <= 0.0))
        {
            return (-2.0 - 3.0 * x);
        }
        else if ((x > 0.0) && (x <= 1.0))
        {
            return (-2.0 + 3.0 * x);
        }
        else if ((x > 1.0) && (x <= 2.0))
        {
            return (2.0 - x);
        }
        else if (x > 2.0)
        {
            return 0;
        }
    }
    return 0; // derivNum outside {0,1,2}; not exercised by any caller, kept as a safe fallback
}

// Solves the tridiagonal system for the natural-cubic-spline control points m_c[1..m_size]
// via Thomas' algorithm (forward elimination then back-substitution), then
// extrapolates the two extra "phantom" control points m_c[0] and m_c[m_size+1]
// needed to evaluate the spline at the very first/last segment.
void CubicSplineInterpolator::solveControlPoints()
{
    double m;
    for (size_t i = 1; i < m_size + 1; i++)
    {
        m_C[i] = 4.0 / 6.0;
        m_A[i] = 1.0 / 6.0;
        m_B[i] = 1.0 / 6.0;
    }
    m_A[1] = 0.0;
    m_B[1] = 0.0;
    m_C[1] = 1.0;
    m_A[m_size] = 0;
    m_B[m_size] = 0;
    m_C[m_size] = m_C[1];
    for (size_t i = 2; i < m_size + 1; i++)
    {
        m = m_A[i] / m_C[i - 1];
        m_C[i] = m_C[i] - m * m_B[i - 1];
        m_x[i] = m_x[i] - m * m_x[i - 1];
    }
    m_c[m_size] = m_x[m_size] / m_C[m_size];
    for (size_t i = m_size - 1; i >= 1; i--)
    {
        m_c[i] = (m_x[i] - m_B[i] * m_c[i + 1]) / m_C[i];
    }
    m_c[0] = 2.0 * m_c[1] - m_c[2];
    m_c[m_size + 1] = 2.0 * m_c[m_size] - m_c[m_size - 1];
}

// Iterative (damped Jacobi/Gauss-Seidel-like) alternative to solveControlPoints(); currently unused
// but kept available since it can be re-enabled from solveControlPoints() for debugging.
void CubicSplineInterpolator::solveControlPointsIterative(size_t itn)
{
    for (size_t n = 0; n < itn; n++)
    {
        m_c[0] = m_c[3] + 3.0 * m_c[1] - 3.0 * m_c[2];
        for (size_t i = 1; i < m_size + 1; i++)
        {
            m_c[i] = m_c[i] * 0.99 + 0.01 * (m_f[i] - m_c[i - 1] * m_A[i] - m_c[i + 1] * m_B[i]) / m_C[i];
        }
        m_c[m_size + 1] = m_c[m_size - 2] - 3.0 * m_c[m_size - 1] + 3.0 * m_c[m_size];

        size_t i = m_size;
        m_c[i] = m_c[i] * 0.99 + 0.01 * (m_f[i] - m_c[i - 1] * m_A[i] - m_c[i + 1] * m_B[i]) / m_C[i];
    }
}

double CubicSplineInterpolator::positionAt(double n)
{
    return evaluateSpline(n, 0);
}

double CubicSplineInterpolator::velocityAt(double n)
{
    return evaluateSpline(n, 1) / 1000 / m_dt;
}

double CubicSplineInterpolator::accelerationAt(double n)
{
    return evaluateSpline(n, 2) / 1000 / m_dt / m_dt;
}

// Evaluates the spline (or its derivative) at fractional index x by summing
// the 4 basis functions with non-zero support at x (indices j-1 .. j+2).
double CubicSplineInterpolator::evaluateSpline(double x, int derivNum )
{
    int j;
    double ff;
    ff = x + 1.0;
    j = (int)ff;
    double sum = 0;
    if (j == (int)m_size)
    {
        j = (int)m_size - 1;
    }

    m_c[0] = m_c[3] + 3.0 * m_c[1] - 3.0 * m_c[2];
    m_c[m_size + 1] = m_c[m_size - 2] - 3.0 * m_c[m_size - 1] + 3.0 * m_c[m_size];
    for (int i = j - 1; i <= j + 2; i++)
    {
        sum += baseSpline(ff - i, derivNum) * m_c[i];
    }
    return sum;
}

// Smoothness-penalized least-squares objective minimized by optimizeByGradientDescent()/optimizeByRandomSearch():
// data-fit error plus lambda * (discrete 3rd-derivative)^2, lambda tuned by the low-pass cutoff `coeff`.
double CubicSplineInterpolator::objectiveValue()
{
    double coeff = 0.3;
    double lambda = (1 / M_PI / coeff) * (1 / M_PI / coeff) * (1 / M_PI / coeff);
    double fRes = 0.0, tmp;
    for (size_t i = 0; i < m_size; i++)
    {
        tmp = evaluateSpline(i) - m_vec[i];
        fRes += tmp * tmp;
    }
    for (size_t i = 1; i < m_size; i++)
    {
        tmp = m_c[i - 1] - 3 * m_c[i] + 3 * m_c[i + 1] - m_c[i + 2];
        fRes += lambda * tmp * tmp;
    }
    return fRes;
}

// Random-search alternative to optimizeByGradientDescent(); currently unused (not called by splineViewer.cpp).
void CubicSplineInterpolator::optimizeByRandomSearch(size_t itn)
{
    double fMin = objectiveValue();
    for (size_t i = 0; i < itn; i++)
    {
        for (size_t j = 0; j <= m_size + 1; j++)
        {
            double cTmp = m_c[j];
            m_c[j] += 0.0001 * (rand() * 1.0 / RAND_MAX - 0.5);
            if (objectiveValue() < fMin)
            {
                fMin = objectiveValue();
            }
            else
            {
                m_c[j] = cTmp;
            }
        }
    }
}

// Finite-difference gradient of objectiveValue() with respect to each control point, used by optimizeByGradientDescent().
void CubicSplineInterpolator::computeObjectiveGradient(double fBefore)
{
    for (size_t j = 1; j <= m_size; j++)
    {
        double dc = 0.1 * m_optimizationStep;
        double cTmp = m_c[j];
        m_c[j] += dc;
        m_cGrad[j] = (objectiveValue() - fBefore) / dc;
        m_c[j] = cTmp;
    }
}

// One gradient-descent step (of `itn` sub-iterations) on objectiveValue(); returns 1 if the
// improvement was still above m_tolerance (i.e. more optimization is worthwhile), else 0.
int CubicSplineInterpolator::optimizeByGradientDescent(size_t itn)
{
    if (fabs(m_dF) > m_tolerance)
    {
        double fBefore = objectiveValue();
        for (size_t i = 0; i < itn; i++)
        {
            computeObjectiveGradient(fBefore);
            for (size_t j = 1; j <= m_size; j++)
            {
                m_c[j] -= 0.001 * m_cGrad[j];
            }
        }
        m_c[0] = m_c[3] + 3.0 * m_c[1] - 3.0 * m_c[2];
        m_c[m_size + 1] = m_c[m_size - 2] - 3.0 * m_c[m_size - 1] + 3.0 * m_c[m_size];
        m_dF = fBefore - objectiveValue();
        return 1;
    }
    return 0;
}