#include "leastsquaressolver.h"
#include "globals.h"

LocalLeastSquaresSolver::LocalLeastSquaresSolver()
{
}

// Fits a Gaussian test function (bench_f) at 39 random points, then checks
// self-test for the weighted-least-squares fit, not used by the interactive viewer.
void LocalLeastSquaresSolver::initializeBenchmark()
{
    SamplePoint3D n;
    FieldDerivatives3D d, d0;
    n.x = 0.0;
    n.y = 0.0;
    n.z = 0.0;
    evaluateBenchmarkDerivatives(n, d0);
    n.f = d0.d[F];

    m_p.clear();
    m_d0.clear();
    m_d.clear();

    m_p.push_back(n);
    m_d0.push_back(d0);
    m_d.push_back(d);

    double sampleRadius = 0.5;
    for (int i = 1; i < 39; i++)
    {
        n.x = sampleRadius * (rand() * 1.0 / RAND_MAX - 0.5);
        n.y = sampleRadius * (rand() * 1.0 / RAND_MAX - 0.5);
        n.z = sampleRadius * (rand() * 1.0 / RAND_MAX - 0.5);
        evaluateBenchmarkDerivatives(n, d0);
        n.f = d0.d[F];
        m_p.push_back(n);
        m_d0.push_back(d0);
        m_d.push_back(d);
    }

    for (int i = 0; i < 39; i++)
    {
        fitDerivatives(m_p[i], m_d[i], 0.5);
    }
}

// Gaussian-eliminates M_ into LU (with partial pivoting, row order tracked in ps[])
// so solveLinearSystem() can be called repeatedly against different right-hand sides.
void LocalLeastSquaresSolver::decomposeLU(void)
{
    int i, j, k, pivotindex = 0;
    double scales[50];
    double normrow, pivot, size, biggest, mult;

    for (i = 0; i < VAR_NUM; i++)
    {
        ps[i] = i; // maps the original row order to the pivoted order
        normrow = 0; // largest magnitude in this row

        for (j = 0; j < VAR_NUM; j++)
        {
            LU[i][j] = M_[i][j];
            if (normrow < fabs(LU[i][j]))
            {
                normrow = fabs(LU[i][j]);
            }
        }
        if (normrow != 0)
        {
            scales[i] = 1.0 / normrow; // for the common row-scaling factors used below
        }
        else
        {
            scales[i] = 0.0;
        }
    }

    // Gaussian elimination with (scaled) partial pivoting
    for (k = 0; k < VAR_NUM - 1; k++)
    {
        biggest = 0;
        for (i = k; i < VAR_NUM; i++)
        {
            size = fabs(LU[ps[i]][k]) * scales[ps[i]];
            if (biggest < size)
            {
                biggest = size;
                pivotindex = i;
            }
        }

        if (biggest == 0)
        {
            pivotindex = 0;
        }

        if (pivotindex != k)
        {
            j = ps[k];
            ps[k] = ps[pivotindex];
            ps[pivotindex] = j;
        }

        pivot = LU[ps[k]][k];

        for (i = k + 1; i < VAR_NUM; i++)
        {
            mult = LU[ps[i]][k] / pivot;
            LU[ps[i]][k] = mult;

            if (mult != 0.0)
            {
                for (j = k + 1; j < VAR_NUM; j++)
                {
                    LU[ps[i]][j] -= mult * LU[ps[k]][j];
                }
            }
        }
    }
}

void LocalLeastSquaresSolver::solveLinearSystem(void)
{
    int i, j;
    double dot;

    for (i = 0; i < VAR_NUM; i++)
    {
        dot = 0;
        for (j = 0; j < i; j++)
        {
            dot += LU[ps[i]][j] * x_m[j];
        }

        x_m[i] = b_m[ps[i]] - dot;
    }

    for (i = VAR_NUM - 1; i >= 0; i--)
    {
        dot = 0.0;

        for (j = i + 1; j < VAR_NUM; j++)
        {
            dot += LU[ps[i]][j] * x_m[j];
        }

        x_m[i] = (x_m[i] - dot) / LU[ps[i]][i];
    }
}

void LocalLeastSquaresSolver::invertNormalMatrix(void)
{
    int i, j;
    decomposeLU();

    for (j = 0; j < VAR_NUM; j++)
    {
        for (i = 0; i < VAR_NUM; i++)
        {
            if (i == j)
            {
                b_m[i] = 1;
            }
            else
            {
                b_m[i] = 0;
            }
        }

        solveLinearSystem();

        for (i = 0; i < VAR_NUM; i++)
        {
            Inv[i][j] = x_m[i];
        }
    }
}

void LocalLeastSquaresSolver::clearMatrix(double m[MAX_EQNS][VAR_NUM])
{
    for (int i = 0; i < 40; i++)
    {
        for (int j = 0; j < 40; j++)
        {
            m[i][j] = 0;
        }
    }
}

void LocalLeastSquaresSolver::clearVector(double v[VAR_NUM])
{
    for (int i = 0; i < 10; i++)
    {
        v[i] = 0;
    }
}

// Full (slow) weighted-least-squares fit: builds the normal equations,
// inverts Mt*W*M explicitly (via invertNormalMatrix()), then multiplies through to get
// the coefficient vector. fitDerivativesFast() below solves the same system
// without forming the explicit inverse, and is what the live code paths use.
void LocalLeastSquaresSolver::fitDerivatives(SamplePoint3D &p, FieldDerivatives3D &res, double delta)
{
    int var_num = VAR_NUM;
    int eq_num = (int)m_p.size();
    // first index is a row number, second index is a column number: M[eq_num][var_num]
    clearMatrix(M_0);
    for (int i = 0; i < eq_num; i++)
    {
        double dxl, dyl, dzl;
        dxl = m_p[i].x - p.x;
        dyl = m_p[i].y - p.y;
        dzl = m_p[i].z - p.z;
        double r2 = dxl * dxl + dyl * dyl + dzl * dzl;
        M_0[i][0] = 1.0;
        M_0[i][1] = dxl;             M_0[i][2] = dyl;       M_0[i][3] = dzl;
        M_0[i][4] = 0.5 * dxl * dxl; M_0[i][5] = dxl * dyl; M_0[i][6] = dxl * dzl;
        M_0[i][7] = 0.5 * dyl * dyl; M_0[i][8] = dyl * dzl; M_0[i][9] = 0.5 * dzl * dzl;

        w[i] = exp(-r2 / (delta * delta));

        b_m[i] = 0.0;
    }

    // Mt*W*M
    clearMatrix(M_);
    for (int i = 0; i < var_num; i++)
    {
        for (int j = 0; j < var_num; j++)
        {
            for (int n = 0; n < eq_num; n++)
            {
                M_[i][j] += M_0[n][i] * M_0[n][j] * w[n];
            }
        }
    }
    invertNormalMatrix(); // (Mt*W*M)^-1 is now in Inv, sized [var_num][var_num]

    for (int i = 0; i < eq_num; i++)
    {
        b_m[i] = m_p[i].f;
    }

    // Mt*W*b
    clearVector(mwb);
    for (int i = 0; i < var_num; i++)
    {
        for (int n = 0; n < eq_num; n++)
        {
            mwb[i] += M_0[n][i] * w[n] * b_m[n];
        }
    }

    clearVector(res.d);
    clearMatrix(MWM);
    for (int i = 0; i < var_num; i++)
    {
        for (int j = 0; j < var_num; j++)
        {
            for (int k = 0; k < var_num; k++)
            {
                MWM[i][j] += M_[i][k] * Inv[k][j];
            }
        }
    }

    for (int i = 0; i < var_num; i++)
    {
        for (int j = 0; j < var_num; j++)
        {
            res.d[i] += Inv[i][j] * mwb[j];
        }
    }
}

// Same weighted least-squares system as fitDerivatives(), solved directly via
void LocalLeastSquaresSolver::fitDerivativesFast(SamplePoint3D &p, FieldDerivatives3D &res, double delta)
{
    int var_num = VAR_NUM;
    int eq_num = (int)m_p.size();

    clearMatrix(M_0);
    for (int i = 0; i < eq_num; i++)
    {
        double dxl, dyl, dzl;
        dxl = m_p[i].x - p.x;
        dyl = m_p[i].y - p.y;
        dzl = m_p[i].z - p.z;
        double r2 = dxl * dxl + dyl * dyl + dzl * dzl;
        M_0[i][0] = 1.0;
        M_0[i][1] = dxl;             M_0[i][2] = dyl;       M_0[i][3] = dzl;
        M_0[i][4] = 0.5 * dxl * dxl; M_0[i][5] = dxl * dyl; M_0[i][6] = dxl * dzl;
        M_0[i][7] = 0.5 * dyl * dyl; M_0[i][8] = dyl * dzl; M_0[i][9] = 0.5 * dzl * dzl;

        w[i] = exp(-r2 / (delta * delta));

        b_m[i] = 0.0;
    }

    clearMatrix(M_);
    for (int i = 0; i < var_num; i++)
    {
        for (int j = 0; j < var_num; j++)
        {
            for (int n = 0; n < eq_num; n++)
            {
                M_[i][j] += M_0[n][i] * M_0[n][j] * w[n];
            }
        }
    }

    for (int i = 0; i < eq_num; i++)
    {
        b_m[i] = m_p[i].f;
    }

    clearVector(mwb);
    for (int i = 0; i < var_num; i++)
    {
        for (int n = 0; n < eq_num; n++)
        {
            mwb[i] += M_0[n][i] * w[n] * b_m[n];
        }
    }

    for (int i = 0; i < var_num; i++)
    {
        b_m[i] = mwb[i];
    }
    decomposeLU();
    solveLinearSystem();

    for (int i = 0; i < var_num; i++)
    {
        res.d[i] = x_m[i];
    }
}

// Same local quadratic fit as fitDerivativesFast(), plus one extra constraint row
// least-squares solve of the discrete Poisson equation at an interior point.
// Only the fitted value (not the full derivative vector) is written back, into p.f.
void LocalLeastSquaresSolver::fitPoissonInterior(SamplePoint3D &p, FieldDerivatives3D &res, double delta)
{
    int var_num = VAR_NUM;
    int eq_num = (int)m_p.size();

    clearMatrix(M_0);
    for (int i = 0; i < eq_num; i++)
    {
        double dxl, dyl, dzl;
        dxl = m_p[i].x - p.x;
        dyl = m_p[i].y - p.y;
        dzl = m_p[i].z - p.z;
        double r2 = dxl * dxl + dyl * dyl + dzl * dzl;
        M_0[i][0] = 1.0;
        M_0[i][1] = dxl;             M_0[i][2] = dyl;       M_0[i][3] = dzl;
        M_0[i][4] = 0.5 * dxl * dxl; M_0[i][5] = dxl * dyl; M_0[i][6] = dxl * dzl;
        M_0[i][7] = 0.5 * dyl * dyl; M_0[i][8] = dyl * dzl; M_0[i][9] = 0.5 * dzl * dzl;

        w[i] = exp(-r2 / (delta * delta));

        b_m[i] = 0.0;
    }
    int i = eq_num; // constraint row for the Laplacian
    eq_num++;
    M_0[i][0] = 0.0;
    M_0[i][1] = 0;             M_0[i][2] = 0;       M_0[i][3] = 0;
    M_0[i][4] = 1;             M_0[i][5] = 0;       M_0[i][6] = 0;
    M_0[i][7] = 1;             M_0[i][8] = 0;       M_0[i][9] = 1;

    w[i] = 1.0;

    b_m[i] = 0.0;

    clearMatrix(M_);
    for (int vi = 0; vi < var_num; vi++)
    {
        for (int j = 0; j < var_num; j++)
        {
            for (int n = 0; n < eq_num; n++)
            {
                M_[vi][j] += M_0[n][vi] * M_0[n][j] * w[n];
            }
        }
    }

    for (int vi = 0; vi < eq_num - 1; vi++)
    {
        b_m[vi] = m_p[vi].f; // from the previous iteration
    }
    b_m[eq_num - 1] = p.rhs; // the Poisson right-hand side

    clearVector(mwb);
    for (int vi = 0; vi < var_num; vi++)
    {
        for (int n = 0; n < eq_num; n++)
        {
            mwb[vi] += M_0[n][vi] * w[n] * b_m[n];
        }
    }

    for (int vi = 0; vi < var_num; vi++)
    {
        b_m[vi] = mwb[vi];
    }

    decomposeLU();
    solveLinearSystem();

    p.f = x_m[0];
}

// Same as fitPoissonInterior(), plus a second constraint row that fixes the
// fitted gradient's component along the wall normal to p.f_bound -- used for
// points on a domain boundary (p.is_boundary >= 0 selects which wall).
void LocalLeastSquaresSolver::fitPoissonBoundary(SamplePoint3D &p, FieldDerivatives3D &res, double delta)
{
    // p is assumed to already be known to be a boundary point (checked by the caller)
    int var_num = VAR_NUM;
    int eq_num = (int)m_p.size();

    clearMatrix(M_0);
    for (int i = 0; i < eq_num; i++)
    {
        double dxl, dyl, dzl;
        dxl = m_p[i].x - p.x;
        dyl = m_p[i].y - p.y;
        dzl = m_p[i].z - p.z;
        double r2 = dxl * dxl + dyl * dyl + dzl * dzl;
        M_0[i][0] = 1.0;
        M_0[i][1] = dxl;             M_0[i][2] = dyl;       M_0[i][3] = dzl;
        M_0[i][4] = 0.5 * dxl * dxl; M_0[i][5] = dxl * dyl; M_0[i][6] = dxl * dzl;
        M_0[i][7] = 0.5 * dyl * dyl; M_0[i][8] = dyl * dzl; M_0[i][9] = 0.5 * dzl * dzl;

        w[i] = exp(-r2 / (delta * delta)); // TODO: make w anisotropic in the wall-normal direction

        b_m[i] = 0.0;
    }
    int i = eq_num; // constraint row for the Laplacian
    eq_num++;
    M_0[i][0] = 0.0;
    M_0[i][1] = 0;             M_0[i][2] = 0;       M_0[i][3] = 0;
    M_0[i][4] = 1;             M_0[i][5] = 0;       M_0[i][6] = 0;
    M_0[i][7] = 1;             M_0[i][8] = 0;       M_0[i][9] = 1;

    w[i] = 1.0;

    b_m[i] = 0.0;

    i = eq_num; // constraint row for the wall-normal gradient
    eq_num++;
    M_0[i][0] = 0.0;
    M_0[i][1] = walls[p.is_boundary].nx; M_0[i][2] = walls[p.is_boundary].ny; M_0[i][3] = walls[p.is_boundary].nz;
    M_0[i][4] = 0.0;                     M_0[i][5] = 0.0;                     M_0[i][6] = 0.0;
    M_0[i][7] = 0.0;                     M_0[i][8] = 0.0;                     M_0[i][9] = 0.0;

    w[i] = 1.0;

    b_m[i] = 0.0;

    clearMatrix(M_);
    for (int vi = 0; vi < var_num; vi++)
    {
        for (int j = 0; j < var_num; j++)
        {
            for (int n = 0; n < eq_num; n++)
            {
                M_[vi][j] += M_0[n][vi] * M_0[n][j] * w[n];
            }
        }
    }

    for (int vi = 0; vi < eq_num - 2; vi++)
    {
        b_m[vi] = m_p[vi].f; // from the previous iteration
    }
    b_m[eq_num - 2] = p.rhs;
    b_m[eq_num - 1] = p.f_bound;

    clearVector(mwb);
    for (int vi = 0; vi < var_num; vi++)
    {
        for (int n = 0; n < eq_num; n++)
        {
            mwb[vi] += M_0[n][vi] * w[n] * b_m[n];
        }
    }

    for (int vi = 0; vi < var_num; vi++)
    {
        b_m[vi] = mwb[vi];
    }
    decomposeLU();
    solveLinearSystem();

    p.f = x_m[0];
}

void LocalLeastSquaresSolver::drawSamplePoints(double sc)
{
    glPointSize(4);
    glBegin(GL_POINTS);
    for (int i = 0; i < (int)m_p.size(); i++)
    {
        glColor3f(sc * 10.0 * m_p[i].f, sc * 10.0 * m_p[i].f, -sc * 10.0 * m_p[i].f);
        glVertex3f(m_p[i].x, m_p[i].y, m_p[i].z);
    }
    glEnd();
}

// ---- benchmark scalar field used only by initializeBenchmark(): f = exp(-r^2) centered at x=0.2 ----

double bench_f(double x, double y, double z)
{
    double r2 = x * x + y * y + z * z;
    return exp(-r2);
}
double bench_dfdx(double x, double y, double z)
{
    double r2 = x * x + y * y + z * z;
    return -2 * x * exp(-r2);
}

double bench_dfdy(double x, double y, double z)
{
    double r2 = x * x + y * y + z * z;
    return -2 * y * exp(-r2);
}

double bench_dfdz(double x, double y, double z)
{
    double r2 = x * x + y * y + z * z;
    return -2 * z * exp(-r2);
}

double bench_d2fdx(double x, double y, double z)
{
    double r2 = x * x + y * y + z * z;
    return (4 * x * x - 2) * exp(-r2);
}

double bench_d2fdy(double x, double y, double z)
{
    double r2 = x * x + y * y + z * z;
    return (4 * y * y - 2) * exp(-r2);
}

double bench_d2fdz(double x, double y, double z)
{
    double r2 = x * x + y * y + z * z;
    return (4 * z * z - 2) * exp(-r2);
}

void LocalLeastSquaresSolver::evaluateBenchmarkDerivatives(SamplePoint3D &p, FieldDerivatives3D &res)
{
    double x, y, z;

    x = p.x - 0.2;
    y = p.y;
    z = p.z;

    res.d[F] = bench_f(x, y, z);
    res.d[FX] = bench_dfdx(x, y, z);
    res.d[FY] = bench_dfdy(x, y, z);
    res.d[FZ] = bench_dfdz(x, y, z);

    res.d[FXX] = bench_d2fdx(x, y, z);
    res.d[FYY] = bench_d2fdy(x, y, z);
    res.d[FZZ] = bench_d2fdz(x, y, z);
}