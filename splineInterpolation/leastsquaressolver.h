#ifndef LEASTSQUARESSOLVER_H
#define LEASTSQUARESSOLVER_H
#include <vector>

// number of fitted coefficients: 1 (value) + 3 (gradient) + 6 (Hessian), see deriv_order below
constexpr int VAR_NUM = 10;

// maximum number of neighbor points + boundary-constraint rows used in one local fit
constexpr int MAX_EQNS = 60;

// one sample point: position (x,y,z), field value f, and (if it lies on a
// domain boundary) the boundary plane index/gradient/Poisson-RHS needed by
struct SamplePoint3D
{
    double x, y, z, f;
    int is_boundary;  // < 0 if this is not a boundary point, otherwise the BoundaryPlane index
    double f_bound;   // prescribed gradient value at the boundary
    double rhs;       // right-hand side of the Poisson equation at this point
};

// the fitted value, gradient, and Hessian of a scalar field at a point (see deriv_order)
struct FieldDerivatives3D
{
    double d[VAR_NUM];
};

// a domain boundary plane: nx*x + ny*y + nz*z = d
struct BoundaryPlane
{
    double nx, ny, nz;
    double d;
};

// Fits a local quadratic model f(x,y,z) ~ value + gradient.(x,y,z) + 0.5*(x,y,z).Hessian.(x,y,z)
// to a weighted point cloud (Gaussian weight by distance) via weighted least squares,
// then reads off the value/gradient/Hessian at the fit's origin. Used to reconstruct
// smooth velocity/acceleration/pressure fields (and their derivatives) from scattered,
// noisy tracked-particle samples.
class LocalLeastSquaresSolver
{
private:
    double M_[MAX_EQNS][VAR_NUM],
           M_0[MAX_EQNS][VAR_NUM],
           MWM[MAX_EQNS][VAR_NUM],
           x_m[MAX_EQNS],
           b_m[MAX_EQNS],
           w[MAX_EQNS],
           mwb[MAX_EQNS],
           LU[MAX_EQNS][VAR_NUM],
           Inv[MAX_EQNS][VAR_NUM];
    int ps[MAX_EQNS];

public:
    // index of each fitted coefficient within FieldDerivatives3D::d / SamplePoint3D
    enum deriv_order
    {
        F,
        FX,
        FY,
        FZ,
        FXX,
        FXY,
        FXZ,
        FYY,
        FYZ,
        FZZ
    };

    LocalLeastSquaresSolver();

    std::vector<SamplePoint3D> m_p;             // point cloud used for the current fit
    std::vector<FieldDerivatives3D> m_d, m_d0;      // fitted / reference (benchmark) derivatives
    std::vector<BoundaryPlane> walls;   // all domain boundary planes

    void initializeBenchmark();
    void decomposeLU(void);
    void solveLinearSystem(void);
    void invertNormalMatrix(void);

    void clearMatrix(double m[MAX_EQNS][VAR_NUM]);
    void clearVector(double v[VAR_NUM]);
    void fitDerivatives(SamplePoint3D &p, FieldDerivatives3D &res, double delta);
    void fitDerivativesFast(SamplePoint3D &p, FieldDerivatives3D &res, double delta);

    void fitPoissonInterior(SamplePoint3D &p, FieldDerivatives3D &res, double delta);
    void fitPoissonBoundary(SamplePoint3D &p, FieldDerivatives3D &res, double delta);

    void drawSamplePoints(double sc);
    void evaluateBenchmarkDerivatives(SamplePoint3D &p, FieldDerivatives3D &res);
};

#endif // LEASTSQUARESSOLVER_H