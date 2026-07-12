#include "globals.h"

int mouseDownX, mouseDownY;
bool isRotating = false;
float rotationX0 = 0.0;
float rotationY0 = 0.0;
float rotationX = rotationX0;
float rotationY = rotationY0;
double mouseX, mouseY;
bool autoRedraw = false;
double displayScale = 500.0;
double viewX = 14.507903;
double viewY = 8.300000;
double viewZ = 24.2;
double forwardX = 0.0;
double forwardY = 0.0;
double forwardZ = 0.0;

int displayMaxI = GRID_NX - 1;
int displayMaxJ = GRID_NY - 1;
int displayMaxK = GRID_NZ - 1;
int currTime = 1;
double zoomScale = 1;

double xmin = 1e10;
double xmax = -1e10;
double ymin = 1e10;
double ymax = -1e10;
double zmin = 1e10;
double zmax = -1e10;

bool drawPoint = true;
bool drawLineSeg = false;
bool drawSpline = true;
bool drawVelocity = false;
bool drawAcceleration = false;
bool drawGridVel = false;
bool drawGridAccel = false;
bool drawGridPressure = false;

double minVel = 1e100;
double maxVel = -1e100;

double minAccel = 1e100;
double maxAccel = -1e100;

double minPressure = 1e100;
double maxPressure = -1e100;

bool rangeCalculated = false;

int isNoOptimized[THREAD_COUNT];