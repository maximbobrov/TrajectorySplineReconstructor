#ifndef TOOLS_H
#define TOOLS_H

#include "globals.h"

double get_time(void);
void get_color(double gval, double minValue, double maxValue);

double ddx_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k, double rhs);
double ddy_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k, double rhs);
double ddz_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k, double rhs);
double d2di_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k);
double d2dj_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k);
double d2dk_a(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k);
double d2di_res(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k);
double d2dj_res(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k);
double d2dk_res(double var_[GRID_NX][GRID_NY][GRID_NZ], int i, int j, int k);
#endif // TOOLS_H