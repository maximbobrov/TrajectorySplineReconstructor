#include "globals.h"
#include "tools.h"
#include "fileloader.h"
#include "splineinterpolator.h"
#include "leastsquaressolver.h"

void optimizeSplineRange(int threadIdx, int startIdx, int endIdx);
void renderFrame(void);
void initializeApplication();
void drawOverlay(void);
void drawSceneAxes(void);

ParticleTrackLoader* trackLoader;

vector<CubicSplineInterpolator*> xSplines;
vector<CubicSplineInterpolator*> ySplines;
vector<CubicSplineInterpolator*> zSplines;

vector<ParticleState> currentParticleStates;

LocalLeastSquaresSolver pressureSolver;

vector<int> particlesGrid[NUM_GRID_CELLS][NUM_GRID_CELLS][NUM_GRID_CELLS];

void drawBitmapText(float x, float y, const char* text, void* font = GLUT_BITMAP_8_BY_13)
{
    glRasterPos2f(x, y);
    for (const char* c = text; *c; ++c)
        glutBitmapCharacter(font, *c);
}

void drawBitmapText3D(float x, float y, float z, const char* text, void* font = GLUT_BITMAP_8_BY_13)
{
    glRasterPos3f(x, y, z);
    for (const char* c = text; *c; ++c)
        glutBitmapCharacter(font, *c);
}

void beginOverlay2D(int width, int height)
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, width, 0, height);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glDisable(GL_DEPTH_TEST);
}

void endOverlay2D()
{
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

void drawColorRamp(float x, float y, float w, float h, double minValue, double maxValue)
{
    glBegin(GL_QUAD_STRIP);
    for (int i = 0; i <= 80; ++i)
    {
        double t = i / 80.0;
        get_color(minValue + (maxValue - minValue) * t, minValue, maxValue);
        glVertex2f(x + w * (float)t, y);
        glVertex2f(x + w * (float)t, y + h);
    }
    glEnd();

    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2f(x, y);
    glVertex2f(x + w, y);
    glVertex2f(x + w, y + h);
    glVertex2f(x, y + h);
    glEnd();
}

void drawSceneAxes(void)
{
    double sx = xmax - xmin;
    double sy = ymax - ymin;
    double sz = zmax - zmin;
    double axisLen = max(max(sx, sy), sz) * 0.22;
    if (axisLen <= 1e-12)
        axisLen = 1.0;

    double x0 = xmin;
    double y0 = ymin;
    double z0 = zmin;

    glLineWidth(3.0f);
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.15f, 0.12f);
    glVertex3f((float)x0, (float)y0, (float)z0);
    glVertex3f((float)(x0 + axisLen), (float)y0, (float)z0);
    glColor3f(0.2f, 1.0f, 0.25f);
    glVertex3f((float)x0, (float)y0, (float)z0);
    glVertex3f((float)x0, (float)(y0 + axisLen), (float)z0);
    glColor3f(0.25f, 0.55f, 1.0f);
    glVertex3f((float)x0, (float)y0, (float)z0);
    glVertex3f((float)x0, (float)y0, (float)(z0 + axisLen));
    glEnd();

    glColor3f(1.0f, 0.15f, 0.12f);
    drawBitmapText3D((float)(x0 + axisLen * 1.06), (float)y0, (float)z0, "X");
    glColor3f(0.2f, 1.0f, 0.25f);
    drawBitmapText3D((float)x0, (float)(y0 + axisLen * 1.06), (float)z0, "Y");
    glColor3f(0.25f, 0.55f, 1.0f);
    drawBitmapText3D((float)x0, (float)y0, (float)(z0 + axisLen * 1.06), "Z");
}

void drawMiniAxis(int width)
{
    float x = 82.0f;
    float y = 82.0f;
    float len = 42.0f;

    glLineWidth(3.0f);
    glBegin(GL_LINES);
    glColor3f(1.0f, 0.15f, 0.12f);
    glVertex2f(x, y);
    glVertex2f(x + len, y);
    glColor3f(0.2f, 1.0f, 0.25f);
    glVertex2f(x, y);
    glVertex2f(x - 24.0f, y + 28.0f);
    glColor3f(0.25f, 0.55f, 1.0f);
    glVertex2f(x, y);
    glVertex2f(x, y + len);
    glEnd();

    glColor3f(1.0f, 0.15f, 0.12f);
    drawBitmapText(x + len + 7.0f, y - 4.0f, "X");
    glColor3f(0.2f, 1.0f, 0.25f);
    drawBitmapText(x - 40.0f, y + 31.0f, "Y");
    glColor3f(0.25f, 0.55f, 1.0f);
    drawBitmapText(x - 4.0f, y + len + 7.0f, "Z");
}

void drawOverlay(void)
{
    int width = glutGet(GLUT_WINDOW_WIDTH);
    int height = glutGet(GLUT_WINDOW_HEIGHT);
    beginOverlay2D(width, height);

    char line[256];
    glColor3f(0.15f, 0.75f, 1.0f);
    drawBitmapText(8.0f, (float)height - 18.0f, "STATUS & DISPLAYED FIELDS");
    glColor3f(0.92f, 0.92f, 0.92f);
    snprintf(line, sizeof(line), "FIELDS  1 points: %s  |  2 raw tracks: %s  |  3 splines: %s  |  5 accel: %s",
             drawPoint ? "on" : "off", drawLineSeg ? "on" : "off", drawSpline ? "on" : "off", drawAcceleration ? "on" : "off");
    drawBitmapText(8.0f, (float)height - 36.0f, line);
    snprintf(line, sizeof(line), "STATE   redraw: %s  |  color gain: %.3g  |  time: %d/%d  |  tracks: %d",
             autoRedraw ? "run" : "paused", displayScale, currTime, trackLoader ? trackLoader->m_timeStepCount : 0,
             trackLoader ? trackLoader->m_splinesCount : 0);
    drawBitmapText(8.0f, (float)height - 54.0f, line);
    snprintf(line, sizeof(line), "CAMERA  pos: %.2f %.2f %.2f  |  rotation: %.1f %.1f  |  zoom: %.3g",
             viewX, viewY, viewZ, rotationX, rotationY, zoomScale);
    drawBitmapText(8.0f, (float)height - 72.0f, line);
    snprintf(line, sizeof(line), "DOMAIN  x %.3g..%.3g  |  y %.3g..%.3g  |  z %.3g..%.3g",
             xmin, xmax, ymin, ymax, zmin, zmax);
    drawBitmapText(8.0f, (float)height - 90.0f, line);

    float rightX = max(520.0f, (float)width - 610.0f);
    glColor3f(1.0f, 0.85f, 0.2f);
    drawBitmapText(rightX, (float)height - 18.0f, "KEYBOARD & MOUSE");
    glColor3f(0.92f, 0.92f, 0.92f);
    drawBitmapText(rightX, (float)height - 36.0f, "DISPLAY  1 points   2 raw tracks   3 splines   [ ] color gain");
    drawBitmapText(rightX, (float)height - 54.0f, "CAMERA   W/S forward-back   A/D strafe   Q/E vertical   mouse drag rotate");
    drawBitmapText(rightX, (float)height - 72.0f, "TIME     9 previous frame   0 next frame   Space continuous redraw");
    drawBitmapText(rightX, (float)height - 90.0f, "TOOLS    V optimize splines   Z calculate errors   X save files");

    double scaleMin = drawAcceleration ? minAccel : minVel;
    double scaleMax = drawAcceleration ? maxAccel : maxVel;
    const char* scaleTitle = drawAcceleration ? "ACCELERATION MAGNITUDE [m^2/s^4]" : "VELOCITY MAGNITUDE [m/s]";
    if (scaleMax <= scaleMin)
    {
        scaleMin = 0.0;
        scaleMax = 1.0;
    }

    float sx = (float)width - 320.0f;
    float sy = 72.0f;
    if (sx < 260.0f)
        sx = 260.0f;
    glColor3f(0.15f, 0.75f, 1.0f);
    drawBitmapText(sx, sy + 52.0f, scaleTitle);
    drawColorRamp(sx, sy + 30.0f, 220.0f, 14.0f, scaleMin, scaleMax);
    glColor3f(0.92f, 0.92f, 0.92f);
    snprintf(line, sizeof(line), "%.3g", scaleMin);
    drawBitmapText(sx, sy + 10.0f, line);
    snprintf(line, sizeof(line), "%.3g", 0.5 * (scaleMin + scaleMax));
    drawBitmapText(sx + 94.0f, sy + 10.0f, line);
    snprintf(line, sizeof(line), "%.3g", scaleMax);
    drawBitmapText(sx + 188.0f, sy + 10.0f, line);

    glColor3f(0.15f, 0.75f, 1.0f);
    drawBitmapText(98.0f, 31.0f, "AXES");
    drawMiniAxis(width);

    endOverlay2D();
}

void renderFrame(void)
{
    double orient_x=0.0;
    double orient_y=0.0;
    double orient_z=5.0;
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glScalef(sc,sc,sc);
    forwardX=cos(-rotationY/20)*cos(-rotationX/20);
    forwardY=cos(-rotationY/20)*sin(-rotationX/20);
    forwardZ=sin(-rotationY/20);
    orient_x=viewX+forwardX;
    orient_y=viewY+forwardY;
    orient_z=viewZ+forwardZ;
    gluLookAt(viewX,viewY,viewZ,orient_x,orient_y,orient_z,0,0,1);

    size_t numVelVec = 20;
    size_t numSplinePoints = 100;
    size_t drawSplineStep = 1;
    size_t drawPointStep = 1;

    if(!rangeCalculated)
    {
        for (size_t s = 0; s < trackLoader->m_splinesCount; s++)
        {
            if(!xSplines[s]->m_moreMinNumber)
                continue;
            for( size_t i=0; i <= numSplinePoints; i++ )
            {
                double vel = xSplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  xSplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints)
                        +    ySplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  ySplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints)
                        +    zSplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  zSplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints);
                if(minVel > sqrt(vel))
                    minVel = sqrt(vel);
                if(maxVel < sqrt(vel))
                    maxVel = sqrt(vel);

                double accel = xSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  xSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints)
                        +      ySplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  ySplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints)
                        +      zSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  zSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints);
                if(minAccel > accel)
                    minAccel = accel;
                if(maxAccel < accel)
                    maxAccel= accel;
            }
        }
        rangeCalculated = true;
    }

    if(drawPoint)
    {
        glPointSize(4);
        glBegin(GL_POINTS);
        for( size_t i = 0; i < currentParticleStates.size(); i++ )
        {
            glColor3f(1,1,1);
            glVertex3f(currentParticleStates[i].x, currentParticleStates[i].y, currentParticleStates[i].z);
        }
        glEnd();
    }

    

    for (size_t s = 0; s < trackLoader->m_splinesCount; s+=drawSplineStep)
    {
        if(!xSplines[s]->m_moreMinNumber)
            continue;
        if(drawLineSeg)
        {
            glLineWidth(2);
            glBegin(GL_LINE_STRIP);
            for( size_t i = 0; i < xSplines[s]->m_size; i+=drawPointStep )
            {
                glColor3f(1.0,1.0,1.0);
                glVertex3f(xSplines[s]->m_vec[i], ySplines[s]->m_vec[i], zSplines[s]->m_vec[i]);
            }
            glEnd();
        }
        if(drawSpline)
        {
            glLineWidth(1);
            glBegin(GL_LINE_STRIP);
            for( size_t i=0; i<=numSplinePoints; i++ )
            {
                double vel = xSplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  xSplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints)
                        +    ySplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  ySplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints)
                        +    zSplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints) *  zSplines[s]->velocityAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints);
                get_color(sqrt(vel), minVel, maxVel);
                glVertex3f(xSplines[s]->positionAt(i*(xSplines[s]->m_size-1)*1.0/numSplinePoints), ySplines[s]->positionAt(i*(ySplines[s]->m_size-1)*1.0/numSplinePoints), zSplines[s]->positionAt(i*(zSplines[s]->m_size-1)*1.0/numSplinePoints));
            }
            glEnd();
        }
        if(drawVelocity)
        {
            
        }
        if(drawAcceleration)
        {
            glLineWidth(2);
            glBegin(GL_LINES);
            for( size_t i=0; i<=numVelVec; i++)
            {
                double accel = xSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec) *  xSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec)
                        +      ySplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec) *  ySplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec)
                        +      zSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec) *  zSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec);
                get_color(accel, minAccel, maxAccel);
                glVertex3f( xSplines[s]->positionAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec), ySplines[s]->positionAt(i*(ySplines[s]->m_size-1)*1.0/numVelVec), zSplines[s]->positionAt(i*(zSplines[s]->m_size-1)*1.0/numVelVec));
                glVertex3f( xSplines[s]->positionAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec)  + displayScale * xSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec)
                            ,ySplines[s]->positionAt(i*(ySplines[s]->m_size-1)*1.0/numVelVec) + displayScale * ySplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec)
                            ,zSplines[s]->positionAt(i*(zSplines[s]->m_size-1)*1.0/numVelVec) + displayScale * zSplines[s]->accelerationAt(i*(xSplines[s]->m_size-1)*1.0/numVelVec));
            }
            glEnd();
        }
    }
    if(drawGridVel)
    {
        for( int k=0; k < min(kNum, GRID_NZ); k+=1)
        {
            for( int i=0; i<min(iNum, GRID_NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j< min(jNum, GRID_NY); j+=1)
                {
                    double vel = trackLoader->m_velField[0][i][j][k] * trackLoader->m_velField[0][i][j][k]
                            +    trackLoader->m_velField[1][i][j][k] * trackLoader->m_velField[1][i][j][k]
                            +    trackLoader->m_velField[2][i][j][k] * trackLoader->m_velField[2][i][j][k];
                    get_color(displayScale * vel, minVel, maxVel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * k * 1.0 / GRID_NZ);

                    double velip = trackLoader->m_velField[0][i+1][j][k] * trackLoader->m_velField[0][i+1][j][k]
                            +      trackLoader->m_velField[1][i+1][j][k] * trackLoader->m_velField[1][i+1][j][k]
                            +      trackLoader->m_velField[2][i+1][j][k] * trackLoader->m_velField[2][i+1][j][k];
                    get_color(displayScale * velip, minVel, maxVel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i+1) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * k * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
        for( int i=0; i<min(iNum, GRID_NX); i+=1)
        {
            for( int k=0; k < min(kNum, GRID_NZ)-1; k+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j<min(jNum, GRID_NY); j+=1)
                {
                    double vel = trackLoader->m_velField[0][i][j][k] * trackLoader->m_velField[0][i][j][k]
                            +    trackLoader->m_velField[1][i][j][k] * trackLoader->m_velField[1][i][j][k]
                            +    trackLoader->m_velField[2][i][j][k] * trackLoader->m_velField[2][i][j][k];
                    get_color(displayScale * vel, minVel, maxVel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);

                    double velkp = trackLoader->m_velField[0][i][j][k+1] * trackLoader->m_velField[0][i][j][k+1]
                            +      trackLoader->m_velField[1][i][j][k+1] * trackLoader->m_velField[1][i][j][k+1]
                            +      trackLoader->m_velField[2][i][j][k+1] * trackLoader->m_velField[2][i][j][k+1];
                    get_color(displayScale * velkp, minVel, maxVel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k+1) * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
        for( int j=0; j<min(jNum, GRID_NY); j+=1)
        {
            for( int i=0; i<min(iNum, GRID_NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int k=0; k < min(kNum, GRID_NZ); k+=1)
                {
                    double vel = trackLoader->m_velField[0][i][j][k] * trackLoader->m_velField[0][i][j][k]
                            +    trackLoader->m_velField[1][i][j][k] * trackLoader->m_velField[1][i][j][k]
                            +    trackLoader->m_velField[2][i][j][k] * trackLoader->m_velField[2][i][j][k];
                    get_color(displayScale * vel, minVel, maxVel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);

                    double velip = trackLoader->m_velField[0][i+1][j][k] * trackLoader->m_velField[0][i+1][j][k]
                            +      trackLoader->m_velField[1][i+1][j][k] * trackLoader->m_velField[1][i+1][j][k]
                            +      trackLoader->m_velField[2][i+1][j][k] * trackLoader->m_velField[2][i+1][j][k];
                    get_color(displayScale * velip, minVel, maxVel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i+1) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
    }

    if(drawGridAccel)
    {
        for( int k=0; k < min(kNum, GRID_NZ); k+=1)
        {
            for( int i=0; i<min(iNum, GRID_NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j< min(jNum, GRID_NY); j+=1)
                {
                    double accel = trackLoader->m_accelField[0][i][j][k] * trackLoader->m_accelField[0][i][j][k]
                            +    trackLoader->m_accelField[1][i][j][k] * trackLoader->m_accelField[1][i][j][k]
                            +    trackLoader->m_accelField[2][i][j][k] * trackLoader->m_accelField[2][i][j][k];
                    get_color(displayScale * accel, minAccel, maxAccel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * k * 1.0 / GRID_NZ);

                    double accelip = trackLoader->m_accelField[0][i+1][j][k] * trackLoader->m_accelField[0][i+1][j][k]
                            +      trackLoader->m_accelField[1][i+1][j][k] * trackLoader->m_accelField[1][i+1][j][k]
                            +      trackLoader->m_accelField[2][i+1][j][k] * trackLoader->m_accelField[2][i+1][j][k];
                    get_color(displayScale * accelip, minAccel, maxAccel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i+1) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * k * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
        for( int i=0; i<min(iNum, GRID_NX); i+=1)
        {
            for( int k=0; k < min(kNum, GRID_NZ)-1; k+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j<min(jNum, GRID_NY); j+=1)
                {
                    double accel = trackLoader->m_accelField[0][i][j][k] * trackLoader->m_accelField[0][i][j][k]
                            +    trackLoader->m_accelField[1][i][j][k] * trackLoader->m_accelField[1][i][j][k]
                            +    trackLoader->m_accelField[2][i][j][k] * trackLoader->m_accelField[2][i][j][k];
                    get_color(displayScale * accel, minAccel, maxAccel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);

                    double accelkp = trackLoader->m_accelField[0][i][j][k+1] * trackLoader->m_accelField[0][i][j][k+1]
                            +      trackLoader->m_accelField[1][i][j][k+1] * trackLoader->m_accelField[1][i][j][k+1]
                            +      trackLoader->m_accelField[2][i][j][k+1] * trackLoader->m_accelField[2][i][j][k+1];
                    get_color(displayScale * accelkp, minAccel, maxAccel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k+1) * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
        for( int j=0; j<min(jNum, GRID_NY); j+=1)
        {
            for( int i=0; i<min(iNum, GRID_NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int k=0; k < min(kNum, GRID_NZ); k+=1)
                {
                    double accel = trackLoader->m_accelField[0][i][j][k] * trackLoader->m_accelField[0][i][j][k]
                            +    trackLoader->m_accelField[1][i][j][k] * trackLoader->m_accelField[1][i][j][k]
                            +    trackLoader->m_accelField[2][i][j][k] * trackLoader->m_accelField[2][i][j][k];
                    get_color(displayScale * accel, minAccel, maxAccel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);

                    double accelip = trackLoader->m_accelField[0][i+1][j][k] * trackLoader->m_accelField[0][i+1][j][k]
                            +      trackLoader->m_accelField[1][i+1][j][k] * trackLoader->m_accelField[1][i+1][j][k]
                            +      trackLoader->m_accelField[2][i+1][j][k] * trackLoader->m_accelField[2][i+1][j][k];
                    get_color(displayScale * accelip, minAccel, maxAccel);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i+1) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
    }

    if(drawGridPressure)
    {
        for( int k=0; k < min(kNum, GRID_NZ); k+=1)
        {
            for( int i=0; i<min(iNum, GRID_NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j< min(jNum, GRID_NY); j+=1)
                {
                    double pressure = trackLoader->m_pressureField[i][j][k];
                    get_color(displayScale * pressure, minPressure, maxPressure);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * k * 1.0 / GRID_NZ);

                    double pressureip = trackLoader->m_pressureField[i+1][j][k];
                    get_color(displayScale * pressureip, minPressure, maxPressure);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i+1) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * k * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
        for( int i=0; i<min(iNum, GRID_NX); i+=1)
        {
            for( int k=0; k < min(kNum, GRID_NZ)-1; k+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int j=0; j<min(jNum, GRID_NY); j+=1)
                {
                    double pressure =  trackLoader->m_pressureField[i][j][k];
                    get_color(displayScale * pressure, minPressure, maxPressure);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);

                    double pressurekp =  trackLoader->m_pressureField[i][j][k+1];
                    get_color(displayScale * pressurekp, minPressure, maxPressure);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k+1) * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
        for( int j=0; j<min(jNum, GRID_NY); j+=1)
        {
            for( int i=0; i<min(iNum, GRID_NX)-1; i+=1)
            {
                glBegin(GL_TRIANGLE_STRIP);
                for( int k=0; k < min(kNum, GRID_NZ); k+=1)
                {
                    double pressure =  trackLoader->m_pressureField[i][j][k];
                    get_color(displayScale * pressure, minPressure, maxPressure);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);

                    double pressureip =  trackLoader->m_pressureField[i+1][j][k];
                    get_color(displayScale * pressureip, minPressure, maxPressure);
                    glVertex3f(DOMAIN_X_MIN + (DOMAIN_X_MAX-DOMAIN_X_MIN) * (i+1) * 1.0 / GRID_NX,   DOMAIN_Y_MIN + (DOMAIN_Y_MAX-DOMAIN_Y_MIN) * (j) * 1.0 / GRID_NY,   DOMAIN_Z_MIN + (DOMAIN_Z_MAX-DOMAIN_Z_MIN) * (k) * 1.0 / GRID_NZ);
                }
                glEnd();
            }
        }
    }

    drawSceneAxes();
    drawOverlay();
    glutSwapBuffers();
    if (autoRedraw==1) glutPostRedisplay();
}

void handleMouseDrag(int x,int y)
{
    if (isRotating==1)
    {
        rotationX=rotationX0+0.1*(x-mouseDownX);
        rotationY=rotationY0+0.1*(y-mouseDownY);
    }
    glutPostRedisplay();
}

void handleMouseButton(int button, int state,int x, int y)
{
    if (state==GLUT_UP)
    {
        isRotating=0;
        rotationX0=rotationX;
        rotationY0=rotationY;
    }
    if (state==GLUT_DOWN)
    {
        isRotating=1;
        mouseDownX=x;
        mouseDownY=y;
    }
    mouseX=(1.0*x)/WINDOW_WIDTH;
    mouseY=(WINDOW_HEIGHT-(1.0*y))/WINDOW_HEIGHT;
    glutPostRedisplay();
}

double loadRuntimeSettings()
{
    char name[5000];
    char name1[5000];
    FILE *file_data=fopen("settings.txt", "r");
    if (file_data == nullptr)
    {
        printf("settings.txt was not found, using built-in defaults\n");
        return 0;
    }
    if (fscanf(file_data,"%4999s %4999s", name, name1) != 2)
    {
        printf("settings.txt has an invalid format, using built-in defaults\n");
    }
    fclose(file_data);
    return 0;
}

void optimizeSplineRange(int threadIdx, int startIdx, int endIdx)
{
    for (int s = startIdx; s < endIdx; s++)
    {
        if(!xSplines[s]->m_moreMinNumber)
            continue;
        isNoOptimized[threadIdx] += xSplines[s]->optimizeByGradientDescent(1);
        isNoOptimized[threadIdx] += ySplines[s]->optimizeByGradientDescent(1);
        isNoOptimized[threadIdx] += zSplines[s]->optimizeByGradientDescent(1);
    }
}

void calculateErrorMetrics()
{
    for (size_t numToSkip = 0; numToSkip < 4; numToSkip++)
    {
        double averPosError = 0.0;
        double averVelError = 0.0;
        double averAccelError = 0.0;
        double meanVel = 0.0;
        int number = 0;
        for (size_t s = 0; s < trackLoader->m_splinesCount; s++)
        {
            if(!xSplines[s]->m_moreMinNumber)
                continue;
            for( size_t i = numToSkip; i < xSplines[s]->m_size - numToSkip; i++ )
            {
                double posX = xSplines[s]->positionAt(i) - trackLoader->m_data[s].x[i];
                double posY = ySplines[s]->positionAt(i) - trackLoader->m_data[s].y[i];
                double posZ = zSplines[s]->positionAt(i) - trackLoader->m_data[s].z[i];
                averPosError += posX * posX + posY * posY + posZ * posZ;

                double vecX = xSplines[s]->velocityAt(i) - trackLoader->m_vel[s].x[i];
                double vecY = ySplines[s]->velocityAt(i) - trackLoader->m_vel[s].y[i];
                double vecZ = zSplines[s]->velocityAt(i) - trackLoader->m_vel[s].z[i];
                averVelError += vecX * vecX + vecY * vecY + vecZ * vecZ;
                meanVel+=sqrt(trackLoader->m_vel[s].x[i]*trackLoader->m_vel[s].x[i] + trackLoader->m_vel[s].y[i]*trackLoader->m_vel[s].y[i] + trackLoader->m_vel[s].z[i]*trackLoader->m_vel[s].z[i]);

                double accelX = xSplines[s]->accelerationAt(i) - trackLoader->m_accel[s].x[i];
                double accelY = ySplines[s]->accelerationAt(i) - trackLoader->m_accel[s].y[i];
                double accelZ = zSplines[s]->accelerationAt(i) - trackLoader->m_accel[s].z[i];
                averAccelError += accelX * accelX + accelY * accelY + accelZ * accelZ;

                number++;
            }
        }
        printf("skip=%d DiffVel = %f DiffAccel = %f DiffPos = %f \n", numToSkip, sqrt(averVelError/number),sqrt(averAccelError/number),sqrt(averPosError/number));
    }
}

void interpolateFieldsToGrid()
{
    printf("Intarpolation started\n");

    double grid_dx = (xmax - xmin) * 1.0 / (NUM_GRID_CELLS);
    double grid_dy = (ymax - ymin) * 1.0 / (NUM_GRID_CELLS);
    double grid_dz = (zmax - zmin) * 1.0 / (NUM_GRID_CELLS);
    double dr = pow(grid_dx*grid_dy*grid_dz,1.0/3) * 1.5;
    for( int i=0; i<GRID_NX; i+=1)
        for( int j=0; j<GRID_NY; j+=1)
            for( int k=0; k<GRID_NZ; k+=1)
            {
                trackLoader->m_numInCell[i][j][k]=0;
                for( int xyz=0; xyz<3; xyz+=1)
                {
                    trackLoader->m_velField[xyz][i][j][k]=0.0;
                    trackLoader->m_accelField[xyz][i][j][k]=0.0;
                }
            }
    for( int i=0; i<GRID_NX; i+=1)
    {
        printf("i=%d\n", i);
        for( int j=0; j<GRID_NY; j+=1)
        {
            for( int k=0; k<GRID_NZ; k+=1)
            {

                trackLoader->m_neighboursInCell[i][j][k].clear();
                int xIdx = int((i*CELL_SIZE_X)/grid_dx);
                int yIdx = int((j*CELL_SIZE_Y)/grid_dy);
                int zIdx = int((k*CELL_SIZE_Z)/grid_dz);
                int im = max(0, xIdx-2);
                int ip = min(NUM_GRID_CELLS - 1, xIdx+2);
                int jm = max(0, yIdx-2);
                int jp = min(NUM_GRID_CELLS - 1, yIdx+2);
                int km = max(0, zIdx-2);
                int kp = min(NUM_GRID_CELLS - 1, zIdx+2);
                double lx = DOMAIN_X_MIN + i*CELL_SIZE_X;
                double ly = DOMAIN_Y_MIN + j*CELL_SIZE_Y;
                double lz = DOMAIN_Z_MIN + k*CELL_SIZE_Z;

                double locdr = dr/2.0;
                while(trackLoader->m_neighboursInCell[i][j][k].size() < 11)
                {
                    locdr*=2.0;
                    trackLoader->m_neighboursInCell[i][j][k].clear();
                    for (int ii = im; ii <= ip ; ii++)
                        for (int jj = jm; jj <=jp; jj++)
                            for (int kk = km; kk <=kp; kk++)
                                for(int n = 0; n < particlesGrid[ii][jj][kk].size(); n++)
                                {
                                    double r2 = (lx - currentParticleStates[particlesGrid[ii][jj][kk].at(n)].x) * (lx - currentParticleStates[particlesGrid[ii][jj][kk].at(n)].x)
                                            +   (ly - currentParticleStates[particlesGrid[ii][jj][kk].at(n)].y) * (ly - currentParticleStates[particlesGrid[ii][jj][kk].at(n)].y)
                                            +   (lz - currentParticleStates[particlesGrid[ii][jj][kk].at(n)].z) * (lz - currentParticleStates[particlesGrid[ii][jj][kk].at(n)].z);
                                    if(r2 <= locdr * locdr)
                                    {
                                        trackLoader->m_neighboursInCell[i][j][k].push_back(particlesGrid[ii][jj][kk].at(n));
                                    }
                                }
                }

                if(trackLoader->m_neighboursInCell[i][j][k].size() > 1)
                    for (int ii = 0; ii < trackLoader->m_neighboursInCell[i][j][k].size()  - 1; ii++)
                        for (int jj = 0; jj < trackLoader->m_neighboursInCell[i][j][k].size()  - ii - 1; jj++)
                        {
                            double r1 = (lx - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj)].x) *   (lx - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj)].x)
                                    +   (ly - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj)].y) *   (ly - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj)].y)
                                    +   (lz - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj)].z) *   (lz - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj)].z);
                            double r2 = (lx - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj+1)].x) * (lx - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj+1)].x)
                                    +   (ly - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj+1)].y) * (ly - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj+1)].y)
                                    +   (lz - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj+1)].z) * (lz - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(jj+1)].z);
                            if (r1 > r2)
                                swap(trackLoader->m_neighboursInCell[i][j][k].at(jj), trackLoader->m_neighboursInCell[i][j][k].at(jj+1));
                        }
            }
        }
    }

    SamplePoint3D n;
    FieldDerivatives3D derivsAx, derivsAy, derivsAz;
    LocalLeastSquaresSolver lssAx, lssAy, lssAz, pressureSolver;

    for( int i=0; i<GRID_NX; i+=1)
    {
        printf("i=%d\n", i);
        for( int j=0; j<GRID_NY; j+=1)
        {
            for( int k=0; k<GRID_NZ; k+=1)
            {
                double lx = DOMAIN_X_MIN + i*CELL_SIZE_X;
                double ly = DOMAIN_Y_MIN + j*CELL_SIZE_Y;
                double lz = DOMAIN_Z_MIN + k*CELL_SIZE_Z;

                if(trackLoader->m_neighboursInCell[i][j][k].size() < 11)
                {
                    printf("lesThenEleven\n");
                    continue;
                }
                if(trackLoader->m_neighboursInCell[i][j][k].size() > 50)
                    trackLoader->m_neighboursInCell[i][j][k].resize(50);

                lssAx.m_p.clear();
                lssAy.m_p.clear();
                lssAz.m_p.clear();
                double rad = 0.5 * ((lx - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(trackLoader->m_neighboursInCell[i][j][k].size() - 1)].x) * (lx - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(trackLoader->m_neighboursInCell[i][j][k].size() - 1)].x)
                        +           (ly - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(trackLoader->m_neighboursInCell[i][j][k].size() - 1)].y) * (ly - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(trackLoader->m_neighboursInCell[i][j][k].size() - 1)].y)
                        +           (lz - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(trackLoader->m_neighboursInCell[i][j][k].size() - 1)].z) * (lz - currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(trackLoader->m_neighboursInCell[i][j][k].size() - 1)].z));
                for (int ii = 0; ii< trackLoader->m_neighboursInCell[i][j][k].size(); ii++)
                {
                    n.x=currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(ii)].x;
                    n.y=currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(ii)].y;
                    n.z=currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(ii)].z;

                    n.f = currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(ii)].ax;
                    lssAx.m_p.push_back(n);
                    n.f = currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(ii)].ay;
                    lssAy.m_p.push_back(n);
                    n.f = currentParticleStates[trackLoader->m_neighboursInCell[i][j][k].at(ii)].az;
                    lssAz.m_p.push_back(n);
                }
                n.x = lx;
                n.y = ly;
                n.z = lz;

                lssAx.fitDerivativesFast(n, derivsAx,  rad);
                lssAy.fitDerivativesFast(n, derivsAy,  rad);
                lssAz.fitDerivativesFast(n, derivsAz,  rad);

                trackLoader->m_accelField[0][i][j][k] = derivsAx.d[LocalLeastSquaresSolver::F];
                trackLoader->m_accelField[1][i][j][k] = derivsAy.d[LocalLeastSquaresSolver::F];
                trackLoader->m_accelField[2][i][j][k] = derivsAz.d[LocalLeastSquaresSolver::F];

                maxAccel = max(maxAccel, sqrt(trackLoader->m_accelField[0][i][j][k]*trackLoader->m_accelField[0][i][j][k]
                        +  trackLoader->m_accelField[1][i][j][k] * trackLoader->m_accelField[1][i][j][k]
                        +  trackLoader->m_accelField[2][i][j][k] * trackLoader->m_accelField[2][i][j][k]));

                minAccel = min(minAccel, sqrt(trackLoader->m_accelField[0][i][j][k]*trackLoader->m_accelField[0][i][j][k]
                        +  trackLoader->m_accelField[1][i][j][k] * trackLoader->m_accelField[1][i][j][k]
                        +  trackLoader->m_accelField[2][i][j][k] * trackLoader->m_accelField[2][i][j][k]));

                trackLoader->m_divAccelField[i][j][k] = derivsAx.d[LocalLeastSquaresSolver::FX]
                        +  derivsAy.d[LocalLeastSquaresSolver::FY]
                        +  derivsAz.d[LocalLeastSquaresSolver::FZ];
            }
        }
    }

    
    printf("Intarpolation finished\n");
}

void smoothGridFields(int iter)
{
    printf("Average started\n");

    for( size_t tt=0; tt < iter; tt++ )
    {
        for( int xyz=0; xyz<3; xyz+=1)
            for( int i=1; i<GRID_NX-1; i+=1)
                for( int j=1; j<GRID_NY-1; j+=1)
                    for( int k = 1; k < GRID_NZ-1; k+=1)
                    {
                        if(trackLoader->m_velFieldCurr[xyz][i][j][k] > 1e-7)
                            trackLoader->m_velField[xyz][i][j][k] = trackLoader->m_velFieldCurr[xyz][i][j][k];
                        if(trackLoader->m_accelFieldCurr[xyz][i][j][k] > 1e-7)
                            trackLoader->m_accelField[xyz][i][j][k] = trackLoader->m_accelFieldCurr[xyz][i][j][k];
                    }

        for( int xyz=0; xyz<3; xyz+=1)
            for( int i=1; i<GRID_NX-1; i+=1)
                for( int j=1; j<GRID_NY-1; j+=1)
                    for( int k = 1; k < GRID_NZ-1; k+=1)
                    {
                        trackLoader->m_velFieldFiltered[xyz][i][j][k] = (trackLoader->m_velField[xyz][i][j][k]
                                                                        + trackLoader->m_velField[xyz][i+1][j][k]
                                + trackLoader->m_velField[xyz][i-1][j][k]
                                + trackLoader->m_velField[xyz][i][j+1][k]
                                + trackLoader->m_velField[xyz][i][j-1][k]
                                + trackLoader->m_velField[xyz][i][j][k+1]
                                + trackLoader->m_velField[xyz][i][j][k-1])/7.0;
                        trackLoader->m_accelFieldFiltered[xyz][i][j][k] = (trackLoader->m_accelField[xyz][i][j][k]
                                                                          + trackLoader->m_accelField[xyz][i+1][j][k]
                                + trackLoader->m_accelField[xyz][i-1][j][k]
                                + trackLoader->m_accelField[xyz][i][j+1][k]
                                + trackLoader->m_accelField[xyz][i][j-1][k]
                                + trackLoader->m_accelField[xyz][i][j][k+1]
                                + trackLoader->m_accelField[xyz][i][j][k-1])/7.0;
                    }

        for( int xyz=0; xyz<3; xyz+=1)
            for( int i=1; i<GRID_NX-1; i+=1)
                for( int j=1; j<GRID_NY-1; j+=1)
                    for( int k = 1; k < GRID_NZ-1; k+=1)
                    {
                        trackLoader->m_velField[xyz][i][j][k] = (trackLoader->m_velFieldFiltered[xyz][i][j][k]
                                                                + trackLoader->m_velFieldFiltered[xyz][i+1][j][k]
                                + trackLoader->m_velFieldFiltered[xyz][i-1][j][k]
                                + trackLoader->m_velFieldFiltered[xyz][i][j+1][k]
                                + trackLoader->m_velFieldFiltered[xyz][i][j-1][k]
                                + trackLoader->m_velFieldFiltered[xyz][i][j][k+1]
                                + trackLoader->m_velFieldFiltered[xyz][i][j][k-1])/7.0;
                        trackLoader->m_accelField[xyz][i][j][k] = (trackLoader->m_accelFieldFiltered[xyz][i][j][k]
                                                                  + trackLoader->m_accelFieldFiltered[xyz][i+1][j][k]
                                + trackLoader->m_accelFieldFiltered[xyz][i-1][j][k]
                                + trackLoader->m_accelFieldFiltered[xyz][i][j+1][k]
                                + trackLoader->m_accelFieldFiltered[xyz][i][j-1][k]
                                + trackLoader->m_accelFieldFiltered[xyz][i][j][k+1]
                                + trackLoader->m_accelFieldFiltered[xyz][i][j][k-1])/7.0;
                    }
    }
    printf("Average finished\n");
}

void updateCurrentParticleStates()
{
    currentParticleStates.clear();
    for (size_t s = 0; s < trackLoader->m_splinesCount; s++)
    {
        if(!xSplines[s]->m_moreMinNumber)
            continue;
        for( size_t i = 0; i < xSplines[s]->m_size; i++ )
        {
            if(currTime == xSplines[s]->m_time[i])
            {
                currentParticleStates.push_back(ParticleState(xSplines[s]->positionAt(i), ySplines[s]->positionAt(i), zSplines[s]->positionAt(i),
                                                 xSplines[s]->velocityAt(i), ySplines[s]->velocityAt(i), zSplines[s]->velocityAt(i),
                                                 xSplines[s]->accelerationAt(i), ySplines[s]->accelerationAt(i), zSplines[s]->accelerationAt(i)));
                currentParticleStates[currentParticleStates.size() - 1].p = rand() * 0.1/ RAND_MAX;
            }
        }
    }
}

void findParticleNeighbors()
{
    double grid_dx = (DOMAIN_X_MAX - DOMAIN_X_MIN) * 1.0 / (NUM_GRID_CELLS);
    double grid_dy = (DOMAIN_Y_MAX - DOMAIN_Y_MIN) * 1.0 / (NUM_GRID_CELLS);
    double grid_dz = (DOMAIN_Z_MAX - DOMAIN_Z_MIN) * 1.0 / (NUM_GRID_CELLS);
    double dr = pow(grid_dx*grid_dy*grid_dz,1.0/3)*2.0;
    for (size_t i = 0; i < NUM_GRID_CELLS; i++)
    {
        for (size_t j = 0; j < NUM_GRID_CELLS; j++)
        {
            for (size_t k = 0; k < NUM_GRID_CELLS; k++)
            {
                particlesGrid[i][j][k].clear();
            }
        }
    }
    for (int i = 0; i < currentParticleStates.size(); i++)
    {
        int xIdx = int((currentParticleStates[i].x - DOMAIN_X_MIN)/grid_dx);
        int yIdx = int((currentParticleStates[i].y - DOMAIN_Y_MIN)/grid_dy);
        int zIdx = int((currentParticleStates[i].z - DOMAIN_Z_MIN)/grid_dz);

        if(xIdx> NUM_GRID_CELLS-1)
            xIdx =NUM_GRID_CELLS-1;

        if(yIdx> NUM_GRID_CELLS-1)
            yIdx =NUM_GRID_CELLS-1;

        if(zIdx> NUM_GRID_CELLS-1)
            zIdx =NUM_GRID_CELLS-1;

        particlesGrid[xIdx][yIdx][zIdx].push_back(i);
    }

    //find boundary by x
    for (int i = 0; i < 2 ; i++)
        for (int j = 0; j < NUM_GRID_CELLS ; j++)
            for (int k = 0; k < NUM_GRID_CELLS; k++)
            {
                int iIdx = i == 0 ? 0 : NUM_GRID_CELLS - 1;
                int size = particlesGrid[iIdx][j][k].size();
                if(size < 2)
                    continue;
                double length = pow (((grid_dx * grid_dy * grid_dz) * 1.0 / size), 1.0/ 3.0);
                int targetNumBound = int(size * length / grid_dx);

                for (int ii = 0; ii < size - 1; ii++)
                    for (int jj = 0; jj < size - ii - 1; jj++)
                        if (currentParticleStates[particlesGrid[iIdx][j][k].at(jj)].x > currentParticleStates[particlesGrid[iIdx][j][k].at(jj+1)].x)
                            swap(currentParticleStates[particlesGrid[iIdx][j][k].at(jj)], currentParticleStates[particlesGrid[iIdx][j][k].at(jj+1)]);

                if(iIdx == 0)
                    for (int n = 0; n < targetNumBound ; n++)
                    {
                        currentParticleStates[particlesGrid[iIdx][j][k].at(n)].isBound = true;
                        currentParticleStates[particlesGrid[iIdx][j][k].at(n)].boundType = 0;
                    }
                else
                    for (int n = size - 1; n > size - targetNumBound ; n--)
                    {
                        currentParticleStates[particlesGrid[iIdx][j][k].at(n)].isBound = true;
                        currentParticleStates[particlesGrid[iIdx][j][k].at(n)].boundType = 1;
                    }
            }

    //find boundary by y
    for (int i = 0; i < NUM_GRID_CELLS ; i++)
        for (int j = 0; j < 2 ; j++)
            for (int k = 0; k < NUM_GRID_CELLS; k++)
            {
                int jIdx = j == 0 ? 0 : NUM_GRID_CELLS - 1;
                int size = particlesGrid[i][jIdx][k].size();
                if(size < 2)
                    continue;
                double length = pow (((grid_dx * grid_dy * grid_dz) * 1.0 / size), 1.0/ 3.0);
                int targetNumBound = int(size * length / grid_dy);

                for (int ii = 0; ii < size - 1; ii++)
                    for (int jj = 0; jj < size - ii - 1; jj++)
                        if (currentParticleStates[particlesGrid[i][jIdx][k].at(jj)].y > currentParticleStates[particlesGrid[i][jIdx][k].at(jj+1)].y)
                            swap(currentParticleStates[particlesGrid[i][jIdx][k].at(jj)], currentParticleStates[particlesGrid[i][jIdx][k].at(jj+1)]);

                if(jIdx == 0)
                    for (int n = 0; n < targetNumBound ; n++)
                    {
                        currentParticleStates[particlesGrid[i][jIdx][k].at(n)].isBound = true;
                        currentParticleStates[particlesGrid[i][jIdx][k].at(n)].boundType = 2;
                    }
                else
                    for (int n = size - 1; n > size - targetNumBound ; n--)
                    {
                        currentParticleStates[particlesGrid[i][jIdx][k].at(n)].isBound = true;
                        currentParticleStates[particlesGrid[i][jIdx][k].at(n)].boundType = 3;
                    }
            }

    //find boundary by z
    for (int i = 0; i < NUM_GRID_CELLS ; i++)
        for (int j = 0; j <  NUM_GRID_CELLS ; j++)
            for (int k = 0; k < 2; k++)
            {
                int kIdx = k == 0 ? 0 : NUM_GRID_CELLS - 1;
                int size = particlesGrid[i][j][kIdx].size();
                if(size < 3)
                    continue;
                double length = grid_dz;//pow (((grid_dx * grid_dy * grid_dz) * 1.0 / size), 1.0/ 3.0);
                int targetNumBound = int(size * length / grid_dz);

                for (int ii = 0; ii < size - 1; ii++)
                    for (int jj = 0; jj < size - ii - 1; jj++)
                        if (currentParticleStates[particlesGrid[i][j][kIdx].at(jj)].z > currentParticleStates[particlesGrid[i][j][kIdx].at(jj+1)].z)
                            swap(currentParticleStates[particlesGrid[i][j][kIdx].at(jj)], currentParticleStates[particlesGrid[i][j][kIdx].at(jj+1)]);

                if(kIdx == 0)
                    for (int n = 0; n < targetNumBound ; n++)
                    {
                        currentParticleStates[particlesGrid[i][j][kIdx].at(n)].isBound = true;
                        currentParticleStates[particlesGrid[i][j][kIdx].at(n)].boundType = 4;
                    }
                else
                    for (int n = size - 1; n > size - targetNumBound ; n--)
                    {
                        currentParticleStates[particlesGrid[i][j][kIdx].at(n)].isBound = true;
                        currentParticleStates[particlesGrid[i][j][kIdx].at(n)].boundType = 5;
                    }
            }

    for (int s = 0; s < currentParticleStates.size(); s++)
    {

        currentParticleStates[s].neighbours.clear();
        int xIdx = int((currentParticleStates[s].x - xmin)/grid_dx);
        int yIdx = int((currentParticleStates[s].y - ymin)/grid_dy);
        int zIdx = int((currentParticleStates[s].z - zmin)/grid_dz);
        int im = max(0, xIdx-2);
        int ip = min(NUM_GRID_CELLS - 1, xIdx+2);
        int jm = max(0, yIdx-2);
        int jp = min(NUM_GRID_CELLS - 1, yIdx+2);
        int km = max(0, zIdx-2);
        int kp = min(NUM_GRID_CELLS - 1, zIdx+2);

        double locdr = dr/2.0;
        while(currentParticleStates[s].neighbours.size() < 11)
        {
            locdr*=2.0;
            currentParticleStates[s].neighbours.clear();
            for (int i = im; i <= ip ; i++)
                for (int j = jm; j <=jp; j++)
                    for (int k = km; k <=kp; k++)
                        for(int n = 0; n < particlesGrid[i][j][k].size(); n++)
                        {
                            if(s == particlesGrid[i][j][k].at(n))
                                continue;
                            double r2 = (currentParticleStates[s].x - currentParticleStates[particlesGrid[i][j][k].at(n)].x) * (currentParticleStates[s].x - currentParticleStates[particlesGrid[i][j][k].at(n)].x)
                                    +   (currentParticleStates[s].y - currentParticleStates[particlesGrid[i][j][k].at(n)].y) * (currentParticleStates[s].y - currentParticleStates[particlesGrid[i][j][k].at(n)].y)
                                    +   (currentParticleStates[s].z - currentParticleStates[particlesGrid[i][j][k].at(n)].z) * (currentParticleStates[s].z - currentParticleStates[particlesGrid[i][j][k].at(n)].z);
                            if(r2 < locdr * locdr)
                            {
                                currentParticleStates[s].neighbours.push_back(particlesGrid[i][j][k].at(n));
                            }
                        }
        }
        if(currentParticleStates[s].neighbours.size() > 1)
            for (int i = 0; i < currentParticleStates[s].neighbours.size() - 1; i++)
                for (int j = 0; j < currentParticleStates[s].neighbours.size() - i - 1; j++)
                {
                    double r1 = (currentParticleStates[s].x - currentParticleStates[currentParticleStates[s].neighbours.at(j)].x) * (currentParticleStates[s].x - currentParticleStates[currentParticleStates[s].neighbours.at(j)].x)
                            +   (currentParticleStates[s].y - currentParticleStates[currentParticleStates[s].neighbours.at(j)].y) * (currentParticleStates[s].y - currentParticleStates[currentParticleStates[s].neighbours.at(j)].y)
                            +   (currentParticleStates[s].z - currentParticleStates[currentParticleStates[s].neighbours.at(j)].z) * (currentParticleStates[s].z - currentParticleStates[currentParticleStates[s].neighbours.at(j)].z);
                    double r2 = (currentParticleStates[s].x - currentParticleStates[currentParticleStates[s].neighbours.at(j+1)].x) * (currentParticleStates[s].x - currentParticleStates[currentParticleStates[s].neighbours.at(j+1)].x)
                            +   (currentParticleStates[s].y - currentParticleStates[currentParticleStates[s].neighbours.at(j+1)].y) * (currentParticleStates[s].y - currentParticleStates[currentParticleStates[s].neighbours.at(j+1)].y)
                            +   (currentParticleStates[s].z - currentParticleStates[currentParticleStates[s].neighbours.at(j+1)].z) * (currentParticleStates[s].z - currentParticleStates[currentParticleStates[s].neighbours.at(j+1)].z);
                    if (r1 > r2)
                        swap(currentParticleStates[s].neighbours.at(j), currentParticleStates[s].neighbours.at(j+1));
                }
    }
}

void saveComputedTracks()
{

    

    
    printf("Saving particles in time %d/%d", currTime, trackLoader->m_timeStepCount);
    char filename[64];
    sprintf(filename, "particles%d_%d.txt", currTime, trackLoader->m_timeStepCount);
    FILE *file_data=fopen(filename,"w");
    for (int s = 0; s < trackLoader->m_splinesCount; s++)
    {
        for( int i = 3; i < xSplines[s]->m_size-3; i++ )
        {
            if(currTime == xSplines[s]->m_time[i])
            {
                fprintf(file_data,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", trackLoader->m_data[s].x[i], trackLoader->m_data[s].y[i], trackLoader->m_data[s].z[i]
                        , xSplines[s]->positionAt(i), ySplines[s]->positionAt(i), zSplines[s]->positionAt(i)
                        , xSplines[s]->velocityAt(i), ySplines[s]->velocityAt(i), zSplines[s]->velocityAt(i)
                        , xSplines[s]->accelerationAt(i), ySplines[s]->accelerationAt(i), zSplines[s]->accelerationAt(i));
            }
        }
    }
    fclose(file_data);
}

void computePressureRightHandSide()
{
    int i, j ,k;

    for (k=1;k<GRID_NZ-1;k++)
    {
        for (j=1;j<GRID_NY-1;j++)
        {
            for (i=1;i<GRID_NX-1;i++)
            {
                trackLoader->m_RHS[i][j][k] = trackLoader->m_divAccelField[i][j][k];
            }
        }
    }
    for (j=1;j<GRID_NY-1;j++)
    {
        for (k=1;k<GRID_NZ-1;k++)
        {
            trackLoader->m_RHS[0][j][k]=trackLoader->m_accelField[0][0][j][k];
            trackLoader->m_RHS[GRID_NX-1][j][k]=trackLoader->m_accelField[0][GRID_NX-1][j][k];
        }
    }

    for (i=1;i<GRID_NX-1;i++)
    {
        for (k=1;k<GRID_NZ-1;k++)
        {
            trackLoader->m_RHS[i][0][k]=trackLoader->m_accelField[1][i][0][k];
            trackLoader->m_RHS[i][GRID_NY-1][k]=trackLoader->m_accelField[1][i][GRID_NY-1][k];
        }
    }

    for (i=1;i<GRID_NX-1;i++)
    {
        for (j=1;j<GRID_NY-1;j++)
        {
            trackLoader->m_RHS[i][j][0]=trackLoader->m_accelField[2][i][j][0];
            trackLoader->m_RHS[i][j][GRID_NZ-1]=trackLoader->m_accelField[2][i][j][GRID_NZ-1];
        }
    }
}

void reconstructParticlePressure()
{
    SamplePoint3D n;
    FieldDerivatives3D derivsAx, derivsAy, derivsAz;
    LocalLeastSquaresSolver lssAx, lssAy, lssAz;
    for (int s = 0; s < currentParticleStates.size(); s++)
    {
        printf("s=%d\n", s);
        if(currentParticleStates[s].neighbours.size() < 11)
            continue;
        if(currentParticleStates[s].neighbours.size() > 50)
            currentParticleStates[s].neighbours.resize(50);

        currentParticleStates[s].hasP = true;
        lssAx.m_p.clear();
        lssAy.m_p.clear();
        lssAz.m_p.clear();
        double rad = 0.5 * ((currentParticleStates[s].x - currentParticleStates[currentParticleStates[s].neighbours.at(currentParticleStates[s].neighbours.size() - 1)].x) * (currentParticleStates[s].x - currentParticleStates[currentParticleStates[s].neighbours.at(currentParticleStates[s].neighbours.size() - 1)].x)
                +           (currentParticleStates[s].y - currentParticleStates[currentParticleStates[s].neighbours.at(currentParticleStates[s].neighbours.size() - 1)].y) * (currentParticleStates[s].y - currentParticleStates[currentParticleStates[s].neighbours.at(currentParticleStates[s].neighbours.size() - 1)].y)
                +           (currentParticleStates[s].z - currentParticleStates[currentParticleStates[s].neighbours.at(currentParticleStates[s].neighbours.size() - 1)].z) * (currentParticleStates[s].z - currentParticleStates[currentParticleStates[s].neighbours.at(currentParticleStates[s].neighbours.size() - 1)].z));
        for (int i = 0; i < currentParticleStates[s].neighbours.size(); i++)
        {
            n.x=currentParticleStates[currentParticleStates[s].neighbours[i]].x;
            n.y=currentParticleStates[currentParticleStates[s].neighbours[i]].y;
            n.z=currentParticleStates[currentParticleStates[s].neighbours[i]].z;
            n.f = currentParticleStates[currentParticleStates[s].neighbours[i]].ax;
            lssAx.m_p.push_back(n);
            n.f = currentParticleStates[currentParticleStates[s].neighbours[i]].ay;
            lssAy.m_p.push_back(n);
            n.f = currentParticleStates[currentParticleStates[s].neighbours[i]].az;
            lssAz.m_p.push_back(n);
        }
        n.x = currentParticleStates[s].x;
        n.y = currentParticleStates[s].y;
        n.z = currentParticleStates[s].z;

        lssAx.fitDerivativesFast(n, derivsAx,  rad);
        lssAy.fitDerivativesFast(n, derivsAy,  rad);
        lssAz.fitDerivativesFast(n, derivsAz,  rad);

        currentParticleStates[s].diva = derivsAx.d[LocalLeastSquaresSolver::FX]
                +  derivsAy.d[LocalLeastSquaresSolver::FY]
                +  derivsAz.d[LocalLeastSquaresSolver::FZ];
    }

    BoundaryPlane bp;

    bp.nx=1.0; bp.ny=0.0; bp.nz=0.0; bp.d=0.0; //0 x0
    pressureSolver.walls.push_back(bp);

    bp.nx=1.0; bp.ny=0.0; bp.nz=0.0; bp.d=1.0; //1 x1
    pressureSolver.walls.push_back(bp);

    bp.nx=0.0; bp.ny=1.0; bp.nz=0.0; bp.d=0.0; //2 y0
    pressureSolver.walls.push_back(bp);

    bp.nx=0.0; bp.ny=1.0; bp.nz=0.0; bp.d=1.0; //3 y2
    pressureSolver.walls.push_back(bp);

    bp.nx=0.0; bp.ny=0.0; bp.nz=1.0; bp.d=0.0; //4 z0
    pressureSolver.walls.push_back(bp);

    bp.nx=0.0; bp.ny=0.0; bp.nz=1.0; bp.d=0.1; //5 z0
    pressureSolver.walls.push_back(bp);

    
    printf("Pressure calculating finished\n");
}

void solvePressurePoissonStep()
{
    SamplePoint3D n;
    for (int s = 0; s < currentParticleStates.size(); s++)
    {
        if(currentParticleStates[s].neighbours.size() < 11)
        {
            printf("lesThenEleven\n");
            continue;
        }
        if(currentParticleStates[s].neighbours.size() > 50)
            currentParticleStates[s].neighbours.resize(50);

        pressureSolver.m_p.clear();

        n.x = currentParticleStates[s].x;
        n.y = currentParticleStates[s].y;
        n.z = currentParticleStates[s].z;
        n.f = currentParticleStates[s].p;
        n.rhs=-currentParticleStates[s].diva;

        if(currentParticleStates[s].isBound)
        {
            n.is_boundary=currentParticleStates[s].boundType;
            n.f_bound=0.0;
            if(n.is_boundary == 0 || n.is_boundary == 1)
                n.f_bound = currentParticleStates[s].ax;
            if(n.is_boundary == 2 || n.is_boundary == 3)
                n.f_bound = currentParticleStates[s].ay;
            if(n.is_boundary == 4 || n.is_boundary == 5)
                n.f_bound = currentParticleStates[s].az;
        }
        else
        {
            n.is_boundary=-1;
            n.f_bound=0.0;
        }
        pressureSolver.m_p.push_back(n);

        for (int i = 0; i < currentParticleStates[s].neighbours.size(); i++)
        {
            n.x = currentParticleStates[currentParticleStates[s].neighbours[i]].x;
            n.y = currentParticleStates[currentParticleStates[s].neighbours[i]].y;
            n.z = currentParticleStates[currentParticleStates[s].neighbours[i]].z;
            n.f = currentParticleStates[currentParticleStates[s].neighbours[i]].p;
            n.rhs=currentParticleStates[currentParticleStates[s].neighbours[i]].diva;
            if(currentParticleStates[currentParticleStates[s].neighbours[i]].isBound)
            {
                n.is_boundary=currentParticleStates[currentParticleStates[s].neighbours[i]].boundType;
                n.f_bound=0.0;
                if(n.is_boundary == 0 || n.is_boundary == 1)
                    n.f_bound = currentParticleStates[currentParticleStates[s].neighbours[i]].ax;
                if(n.is_boundary == 2 || n.is_boundary == 3)
                    n.f_bound = currentParticleStates[currentParticleStates[s].neighbours[i]].ay;
                if(n.is_boundary == 4 || n.is_boundary == 5)
                    n.f_bound = currentParticleStates[currentParticleStates[s].neighbours[i]].az;
            }
            else
            {
                n.is_boundary=-1;
                n.f_bound=0.0;
            }
            pressureSolver.m_p.push_back(n);
        }
        if (pressureSolver.m_p[0].is_boundary>=0)
        {
            pressureSolver.fitPoissonBoundary(pressureSolver.m_p[0],pressureSolver.m_d[0],0.125);

        }else
        {
            pressureSolver.fitPoissonInterior(pressureSolver.m_p[0],pressureSolver.m_d[0],0.125);
        }
        currentParticleStates[s].p =  pressureSolver.m_p[0].f;
    }
    printf("done\n");

}

void solvePressurePoissonIterations(int in_) //poisson equation at the cell centers
{
    int i,j,k,nn;
    double a,b,c,d,axb_max;
    a=CELL_SIZE_X*CELL_SIZE_X*2.0*(1.0/(CELL_SIZE_X*CELL_SIZE_X)+ 1.0/(CELL_SIZE_Y*CELL_SIZE_Y)+1.0/(CELL_SIZE_Z*CELL_SIZE_Z));
    b=-1.0*CELL_SIZE_X*CELL_SIZE_X/(CELL_SIZE_X*CELL_SIZE_X);
    c=-1.0*CELL_SIZE_X*CELL_SIZE_X/(CELL_SIZE_Y*CELL_SIZE_Y);
    d=-1.0*CELL_SIZE_X*CELL_SIZE_X/(CELL_SIZE_Z*CELL_SIZE_Z);

    axb_max=0.0;
    for (nn=0;nn<in_;nn++)
    {
        computePressureRightHandSide();
        //neumann bc
        for (j=1;j<GRID_NY-1;j++)
        {
            for (k=1;k<GRID_NZ-1;k++)
            {
                trackLoader->m_pressureField[0][j][k]=ddx_a(trackLoader->m_pressureField,0,j,k,trackLoader->m_RHS[0][j][k]);//0
                trackLoader->m_pressureField[GRID_NX-1][j][k]=ddx_a(trackLoader->m_pressureField,GRID_NX-1,j,k,trackLoader->m_RHS[GRID_NX-1][j][k]);//0
            }
        }

        for (i=1;i<GRID_NX-1;i++)
        {
            for (k=1;k<GRID_NZ-1;k++)
            {
                trackLoader->m_pressureField[i][0][k]=ddy_a(trackLoader->m_pressureField,i,0,k,trackLoader->m_RHS[i][0][k]);
                trackLoader->m_pressureField[i][GRID_NY-1][k]=ddy_a(trackLoader->m_pressureField,i,GRID_NY-1,k,trackLoader->m_RHS[i][GRID_NY-1][k]);
            }
        }

        for (i=1;i<GRID_NX-1;i++)
        {
            for (j=1;j<GRID_NY-1;j++)
            {
                trackLoader->m_pressureField[i][j][0]=ddz_a(trackLoader->m_pressureField,i,j,0,trackLoader->m_RHS[i][j][0]);//0
                trackLoader->m_pressureField[i][j][GRID_NZ-1]=ddz_a(trackLoader->m_pressureField,i,j,GRID_NZ-1,trackLoader->m_RHS[i][j][GRID_NZ-1]);//0
            }
        }

        for (i=1;i<GRID_NX-1;i++)
        {
            for (j=1;j<GRID_NY-1;j++)
            {
                for (k=1;k<GRID_NZ-1;k++)
                {
                    double d2phi_dx_a,d2phi_dy_a,d2phi_dz_a,d2phi_dx_res,d2phi_dy_res,d2phi_dz_res,phi_a,phi_res;

                    d2phi_dx_a=d2di_a(trackLoader->m_pressureField,i,j,k)/(CELL_SIZE_X*CELL_SIZE_X);
                    d2phi_dy_a=d2dj_a(trackLoader->m_pressureField,i,j,k)/(CELL_SIZE_Y*CELL_SIZE_Y);
                    d2phi_dz_a=d2dk_a(trackLoader->m_pressureField,i,j,k)/(CELL_SIZE_Z*CELL_SIZE_Z);

                    d2phi_dx_res=d2di_res(trackLoader->m_pressureField,i,j,k)/(CELL_SIZE_X*CELL_SIZE_X);
                    d2phi_dy_res=d2dj_res(trackLoader->m_pressureField,i,j,k)/(CELL_SIZE_Y*CELL_SIZE_Y);
                    d2phi_dz_res=d2dk_res(trackLoader->m_pressureField,i,j,k)/(CELL_SIZE_Z*CELL_SIZE_Z);

                    phi_a=d2phi_dx_a+d2phi_dy_a+d2phi_dz_a;
                    phi_res=d2phi_dx_res+d2phi_dy_res+d2phi_dz_res;
                    trackLoader->m_pressureField[i][j][k]=-(phi_res+trackLoader->m_RHS[i][j][k])/phi_a;
                }
            }
        }

        double pressure_mean=0.0;
        int nnn=0;

        for (i=1;i<GRID_NX-1;i++)
        {
            for (j=1;j<GRID_NY-1;j++)
            {
                for (k=1;k<GRID_NZ-1;k++)
                {
                    pressure_mean += trackLoader->m_pressureField[i][j][k];
                    nnn++;
                }
            }
        }

        pressure_mean/=nnn;
        for (i=1;i<GRID_NX-1;i++)
        {
            for (j=1;j<GRID_NY-1;j++)
            {
                for (k=1;k<GRID_NZ-1;k++)
                {
                    trackLoader->m_pressureField[i][j][k]-=pressure_mean;
                    if(trackLoader->m_pressureField[i][j][k]< minPressure)
                        minPressure = trackLoader->m_pressureField[i][j][k];
                    if(trackLoader->m_pressureField[i][j][k]> maxPressure)
                        maxPressure = trackLoader->m_pressureField[i][j][k];
                }
            }
        }
    }
}

void handleKeyboard(unsigned char key, int x, int y)
{
    if (key=='1')
    {
        drawPoint = !drawPoint;
    }
    if (key=='2')
    {
        drawLineSeg = !drawLineSeg;
    }
    if (key=='3')
    {
        drawSpline = !drawSpline;
    }
    if (key=='4')
    {
        drawVelocity = !drawVelocity;
    }
    if (key=='5')
    {
        drawAcceleration = !drawAcceleration;
    }
    if (key=='0')
    {
        if(currTime < trackLoader->m_timeStepCount)
        {
            currTime+=1;
            updateCurrentParticleStates();
        }
        printf("Time %d/%d\n", currTime, trackLoader->m_timeStepCount);
    }
    if (key=='9')
    {
        if(currTime > 1)
        {
            currTime-=1;
            updateCurrentParticleStates();
        }
        printf("Time %d/%d\n", currTime, trackLoader->m_timeStepCount);
    }
    if (key=='b')
    {
        smoothGridFields(100);
    }
    if (key=='o')
    {
        iNum++;
        if(iNum>=GRID_NX)
            iNum=GRID_NX-1;
    }
    if (key=='i')
    {
        iNum--;
        if(iNum<0)
            iNum=0;
    }
    if (key=='l')
    {
        jNum++;
        if(jNum>=GRID_NY)
            jNum=GRID_NY-1;
    }
    if (key=='k')
    {
        jNum--;
        if(jNum<0)
            jNum=0;
    }
    if (key=='m')
    {
        for(int i = 0; i<20 ; i++)
            solvePressurePoissonStep();
    }

    if (key=='.')
    {
        kNum++;
        if(kNum>=GRID_NZ)
            kNum=GRID_NZ-1;
    }
    if (key==',')
    {
        kNum--;
        if(kNum<0)
            kNum=0;
    }
    if (key=='[')
    {
        displayScale/=1.2;
    }
    if (key==']')
    {
        displayScale*=1.2;
    }
    if (key=='w')
    {
        viewX+=(forwardX)*1.5;
        viewY+=(forwardY)*1.5;
        viewZ+=(forwardZ)*1.5;
    }
    if (key=='s')
    {
        viewX-=(forwardX)*1.5;
        viewY-=(forwardY)*1.5;
        viewZ-=(forwardZ)*1.5;
    }
    if (key=='q')
    {
        viewZ+=1.1;
    }
    if (key=='e')
    {
        viewZ-=1.1;
    }
    if (key=='a')
    {
        double l2=sqrt(forwardY*forwardY+forwardX*forwardX);
        viewY+=(forwardX)*1.5/l2;
        viewX+=-(forwardY)*1.5/l2;
    }
    if (key=='d')
    {
        double    l2=sqrt(forwardY*forwardY+forwardX*forwardX);
        viewX+=(forwardY)*1.5/l2;
        viewY+=-(forwardX)*1.5/l2;
    }

    if(key == 'z')
    {
        calculateErrorMetrics();
    }
    if(key == 'x')
    {
        saveComputedTracks();
    }
    if(key == 'v')
    {
        int numNoMin = 1;
        int numToStop = 0;
        double time1=get_time();
        while (numNoMin != 0 && numToStop < 10000)
        {
            numToStop++;
            numNoMin = 0;
            std::vector <std::thread> th_vec;
            th_vec.clear();
            for (int i = 0; i < THREAD_COUNT; ++i)
            {
                int numForOneProc = THREAD_COUNT == 1 ? trackLoader->m_splinesCount : (int)(trackLoader->m_splinesCount / (THREAD_COUNT - 1));
                int startIdx = i * numForOneProc;
                int endIdx =  (i + 1) * numForOneProc;
                if(i == (THREAD_COUNT - 1))
                    endIdx = trackLoader->m_splinesCount;
                isNoOptimized[i] = 0;
                th_vec.push_back(std::thread(optimizeSplineRange,i, startIdx, endIdx));
            }
            for (int i = 0; i < THREAD_COUNT; ++i)
            {
                th_vec.at(i).join();
                numNoMin += isNoOptimized[i];
            }
            printf("Number no minimazed splines = %d\n", numNoMin);
        }
        double time2=get_time();
        printf("Optimization time = %e \n", time2-time1);
    }
    if (key==' ')
    {
        autoRedraw=!autoRedraw;
    }

    glutPostRedisplay();
}

void initializeApplication()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(45.0f, WINDOW_WIDTH*1.0/WINDOW_HEIGHT, 0.1f, 250.0f);
    glMatrixMode (GL_MODELVIEW);
    loadRuntimeSettings();
    trackLoader = new ParticleTrackLoader(DATA_DIRECTORY);
    trackLoader->findDataFiles();
    trackLoader->loadParticleTracks();
    printf("%d tracks loaded\n", trackLoader->m_splinesCount);
    xSplines.resize(trackLoader->m_splinesCount);
    ySplines.resize(trackLoader->m_splinesCount);
    zSplines.resize(trackLoader->m_splinesCount);
    int numbers = 0;
    for (size_t i = 0; i < trackLoader->m_splinesCount; i+=1)
    {

        xSplines[numbers] = (new CubicSplineInterpolator(trackLoader->m_data[i].x, trackLoader->m_data[i].t));
        ySplines[numbers] = (new CubicSplineInterpolator(trackLoader->m_data[i].y, trackLoader->m_data[i].t));
        zSplines[numbers] = (new CubicSplineInterpolator(trackLoader->m_data[i].z, trackLoader->m_data[i].t));
        numbers ++;
    }
    trackLoader->m_splinesCount = numbers;
    xSplines.resize(numbers);
    ySplines.resize(numbers);
    zSplines.resize(numbers);
    printf("Splines calculated\n");

    updateCurrentParticleStates();

    viewX=xmin-0.5*(xmax-xmin);
    viewY=0.5*(ymax-ymin);
    viewZ=0.5*(zmax-zmin);
}

int main(int argc, char** argv)
{
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("TrajectorySplineReconstructor");
    glutDisplayFunc(renderFrame);
    glutMotionFunc(handleMouseDrag);
    glutMouseFunc(handleMouseButton);
    glutKeyboardFunc(handleKeyboard);
    initializeApplication();
    glutMainLoop();

    delete trackLoader;
    for (size_t i = 0; i < xSplines.size(); i++) {
        delete xSplines[i];
        delete ySplines[i];
        delete zSplines[i];
    }
    return 0;
}
