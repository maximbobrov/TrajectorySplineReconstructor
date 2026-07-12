#include "tools.h"

// Wall-clock time in seconds (fractional). Only used to measure elapsed
// durations (see the 'v' full-optimization command in splineViewer.cpp), so
// an arbitrary epoch is fine -- replaces the original POSIX gettimeofday()
// (<sys/time.h>, unavailable on MSVC) with std::chrono.
double get_time(void)
{
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// Maps gval (clamped to [minValue, maxValue]) to an RGB color via a 5-stop
// blue -> cyan -> green -> yellow -> red gradient and sets it as the current OpenGL color.
void get_color(double gval, double minValue, double maxValue)
{
    constexpr int NUM_SEGMENTS = 4;
    double val = gval;
    if (val > maxValue)
    {
        val = maxValue;
    }
    if (val < minValue)
    {
        val = minValue;
    }

    struct Xyz
    {
        double x, y, z;
    };

    Xyz colorTable[5];
    colorTable[0] = { 0.0, 0.0, 1.0 };
    colorTable[1] = { 0.0, 1.0, 1.0 };
    colorTable[2] = { 0.0, 1.0, 0.0 };
    colorTable[3] = { 1.0, 1.0, 0.0 };
    colorTable[4] = { 1.0, 0.0, 0.0 };

    int i;
    double alpha;
    if ((maxValue - minValue) > 1e-35)
    {
        alpha = (val - minValue) / (maxValue - minValue) * NUM_SEGMENTS;
        i = (int)(alpha);
        alpha = alpha - i;
    }
    else
    {
        alpha = 0.0;
        i = 2;
    }
    glColor3f(colorTable[i].x * (1 - alpha) + colorTable[i + 1].x * alpha,
              colorTable[i].y * (1 - alpha) + colorTable[i + 1].y * alpha,
              colorTable[i].z * (1 - alpha) + colorTable[i + 1].z * alpha);
}

// The four functions below (ddx_a/ddy_a/ddz_a + d2d*_a) return the coefficient
// of the Neumann (zero-gradient) boundary extrapolation formula var[boundary] = (rhs term)/a
// used by solvePressurePoissonIterations(); the four d2d*_res functions return the corresponding
// interior-neighbor sum. i==0/i==GRID_NX-1 (and the j/k equivalents) are the
// domain boundaries, using a second-order one-sided difference; everything
// else falls back to the plain second-difference stencil.

double ddx_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k, double rhs)
{
    double res;
    if (i == 0)
    {
        res = (-rhs * 2.0 * CELL_SIZE_X + 4.0 * var_[i + 1][j][k] - var_[i + 2][j][k]) / 3.0;
    }
    else
    {
        if (i == GRID_NX - 1)
        {
            res = (rhs * 2.0 * CELL_SIZE_X + 4.0 * var_[i - 1][j][k] - var_[i - 2][j][k]) / 3.0;
        }
        else
        {
            res = 0;
        }
    }
    return res;
}

double ddy_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k, double rhs)
{
    double res;
    if (j == 0)
    {
        res = (-rhs * 2.0 * CELL_SIZE_Y + 4.0 * var_[i][j + 1][k] - var_[i][j + 2][k]) / 3.0;
    }
    else
    {
        if (j == GRID_NY - 1)
        {
            res = (rhs * 2.0 * CELL_SIZE_Y + 4.0 * var_[i][j - 1][k] - var_[i][j - 2][k]) / 3.0;
        }
        else
        {
            res = 0;
        }
    }
    return res;
}

double ddz_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k, double rhs)
{
    double res;
    if (k == 0)
    {
        res = (-rhs * 2.0 * CELL_SIZE_Z + 4.0 * var_[i][j][k + 1] - var_[i][j][k + 2]) / 3.0;
    }
    else
    {
        if (k == GRID_NZ - 1)
        {
            res = (rhs * 2.0 * CELL_SIZE_Z + 4.0 * var_[i][j][k - 1] - var_[i][j][k - 2]) / 3.0;
        }
        else
        {
            res = 0;
        }
    }
    return res;
}

double d2di_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k)
{
    double res;
    if (i == 0)
    {
        res = 1.0;
    }
    else
    {
        if (i == GRID_NX - 1)
        {
            res = 1.0;
        }
        else
        {
            res = -2.0;
        }
    }
    return res;
}

double d2dj_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k)
{
    double res;
    if (j == 0)
    {
        res = 1.0;
    }
    else
    {
        if (j == GRID_NY - 1)
        {
            res = 1.0;
        }
        else
        {
            res = -2.0;
        }
    }
    return res;
}

double d2dk_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k)
{
    double res;
    if (k == 0)
    {
        res = 1.0;
    }
    else
    {
        if (k == GRID_NZ - 1)
        {
            res = 1.0;
        }
        else
        {
            res = -2.0;
        }
    }
    return res;
}

double d2di_res(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k)
{
    double res;
    if (i == 0)
    {
        res = (-2.0 * var_[i + 1][j][k] + var_[i + 2][j][k]);
    }
    else
    {
        if (i == GRID_NX - 1)
        {
            res = (-2.0 * var_[i - 1][j][k] + var_[i - 2][j][k]);
        }
        else
        {
            res = (var_[i + 1][j][k] + var_[i - 1][j][k]);
        }
    }
    return res;
}

double d2dj_res(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k)
{
    double res;
    if (j == 0)
    {
        res = (-2.0 * var_[i][j + 1][k] + var_[i][j + 2][k]);
    }
    else
    {
        if (j == GRID_NY - 1)
        {
            res = (-2.0 * var_[i][j - 1][k] + var_[i][j - 2][k]);
        }
        else
        {
            res = (var_[i][j + 1][k] + var_[i][j - 1][k]);
        }
    }
    return res;
}

double d2dk_res(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k)
{
    double res;
    if (k == 0)
    {
        res = (-2.0 * var_[i][j][k + 1] + var_[i][j][k + 2]);
    }
    else
    {
        if (k == GRID_NZ - 1)
        {
            res = (-2.0 * var_[i][j][k - 1] + var_[i][j][k - 2]);
        }
        else
        {
            res = (var_[i][j][k + 1] + var_[i][j][k - 1]);
        }
    }
    return res;
}