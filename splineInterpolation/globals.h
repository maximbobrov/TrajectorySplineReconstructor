#ifndef GLOBALS_H
#define GLOBALS_H

// Windows.h defines min/max macros that clash with std::min/std::max used
// throughout this project (e.g. tools.cpp, main.cpp) -- NOMINMAX suppresses them.
#define NOMINMAX
// M_PI is not defined by <math.h> on MSVC unless this is defined first.
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include "my_include/gl.h"
#include "my_include/glu.h"
#include "my_include/glut.h"
#include <math.h>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <thread>

// grid resolution of the (r,z,phi)-like Cartesian field grid used by the
// (currently unreachable, see README "Notes") grid-interpolation/pressure-solve code path
constexpr int GRID_NX = 251;
constexpr int GRID_NY = 126;
constexpr int GRID_NZ = 76;

// resolution of the spatial hash grid used to find nearest-neighbour tracked particles
constexpr int NUM_GRID_CELLS = 20;

constexpr int WINDOW_WIDTH = 1600;
constexpr int WINDOW_HEIGHT = 1300;

// physical domain bounds covered by the field grid, meters
constexpr double DOMAIN_X_MIN = 1.75;
constexpr double DOMAIN_X_MAX = 6.75;
constexpr double DOMAIN_Y_MIN = -1.25;
constexpr double DOMAIN_Y_MAX = 1.25;
constexpr double DOMAIN_Z_MIN = 0.04;
constexpr double DOMAIN_Z_MAX = 1.0;

constexpr double CELL_SIZE_X = (DOMAIN_X_MAX - DOMAIN_X_MIN) * 1.0 / GRID_NX;
constexpr double CELL_SIZE_Y = (DOMAIN_Y_MAX - DOMAIN_Y_MIN) * 1.0 / GRID_NY;
constexpr double CELL_SIZE_Z = (DOMAIN_Z_MAX - DOMAIN_Z_MIN) * 1.0 / GRID_NZ;

constexpr int THREAD_COUNT = 16; // number of worker threads used by the 'v' full-optimization command

// folder (relative to the working directory) scanned for input ".dat" track files
constexpr const char* DATA_DIRECTORY = "DA_ppp_0_005/";

using namespace std;

// per-particle sample used for the "current time step" snapshot drawn on screen
struct ParticleState
{
    ParticleState(double ix, double iy, double iz, double iu, double iv, double iw, double iax, double iay, double iaz)
        : x(ix), y(iy), z(iz), u(iu), v(iv), w(iw), ax(iax), ay(iay), az(iaz) {}
    double x, y, z;
    double u, v, w;
    double ax, ay, az, diva;
    bool hasP = false;
    double p;
    vector<int> neighbours;
    bool isBound = false;
    int boundType = -1;
};

// ---- view / camera state (see splineViewer.cpp: onMouseMotion, onMouseButton, onKeyboard) ----
extern int mouseDownX, mouseDownY; // mouse position latched on button-down, for drag deltas
extern bool isRotating;            // true while the left mouse button is held and dragged
extern float rotationX0, rotationY0; // rotation angles latched on button-down
extern float rotationX, rotationY;   // current view rotation angles (updated by dragging)
extern double mouseX, mouseY;        // last mouse position in normalized [0,1] window coordinates
extern bool autoRedraw;              // toggled by Space: keep calling glutPostRedisplay every frame
extern double displayScale;          // color/vector-length multiplier, keys '[' and ']'
extern double viewX, viewY, viewZ;   // camera position
extern double forwardX, forwardY, forwardZ; // camera forward (look) direction, unit vector

extern int displayMaxI; // upper index bound when drawing the field grid (keys 'i'/'o')
extern int displayMaxJ; // upper index bound when drawing the field grid (keys 'k'/'l')
extern int displayMaxK; // upper index bound when drawing the field grid (keys ','/'.')
extern int currTime;    // currently displayed time step (0-based file index), keys '9'/'0'
extern double zoomScale; // scene scale factor applied in onDisplay (glScalef)

#define iNum displayMaxI
#define jNum displayMaxJ
#define kNum displayMaxK
#define sc zoomScale

extern double xmin;
extern double xmax;
extern double ymin;
extern double ymax;
extern double zmin;
extern double zmax;

extern bool drawPoint;         // key '1': draw current-time particle positions
extern bool drawLineSeg;       // key '2': draw raw per-track polylines
extern bool drawSpline;        // key '3': draw fitted splines, colored by speed
extern bool drawVelocity;      // key '4' (see README "Notes": currently has no visible effect)
extern bool drawAcceleration;  // key '5': draw per-track acceleration vectors
extern bool drawGridVel;       // grid velocity-field overlay (see README "Notes": currently unreachable)
extern bool drawGridAccel;     // grid acceleration-field overlay (see README "Notes": currently unreachable)
extern bool drawGridPressure;  // grid pressure-field overlay (see README "Notes": currently unreachable)

extern double minVel;
extern double maxVel;

extern double minAccel;
extern double maxAccel;

extern double minPressure;
extern double maxPressure;

extern bool rangeCalculated; // true once the min/max color-scale ranges have been computed once

extern int isNoOptimized[THREAD_COUNT]; // per-thread "still improving" flag used by the 'v' command

#endif // GLOBALS_H