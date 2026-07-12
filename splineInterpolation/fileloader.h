#ifndef FILELOADER_H
#define FILELOADER_H
#include "globals.h"

// Scans a directory for ".dat" particle-track files (one file per recorded time
// step) and loads them into per-track (x,y,z,t) samples (m_data), keyed by the
// numeric particle id used in the file so samples for the same particle across
// different time-step files are appended to the same track. Also owns the
// (currently unreachable from the UI, see README "Notes") field-grid buffers
// used by the grid velocity/acceleration/pressure interpolation code.
struct ParticleTrackLoader
{
    // one track: position samples (x,y,z) and the time-step index (t) each was recorded at
    struct ParticleTrack
    {
        ParticleTrack();
        ParticleTrack(double ix, double iy, double iz, int it);
        vector<double> x;
        vector<double> y;
        vector<double> z;
        vector<int> t;
    };

    ParticleTrackLoader(string iPath);
    void findDataFiles();
    void loadParticleTracks();

    string m_path;
    int m_numbers[10000000]; // particle id -> index into m_data, or -1 if not seen yet
    vector<ParticleTrack> m_data;

    // field grid used by the grid interpolation / pressure-solve code path (see README "Notes")
    int m_numInCell[GRID_NX][GRID_NY][GRID_NZ];
    vector<int> m_neighboursInCell[GRID_NX][GRID_NY][GRID_NZ];
    double m_velField[3][GRID_NX][GRID_NY][GRID_NZ];
    double m_velFieldCurr[3][GRID_NX][GRID_NY][GRID_NZ];
    double m_velFieldFiltered[3][GRID_NX][GRID_NY][GRID_NZ];
    double m_accelField[3][GRID_NX][GRID_NY][GRID_NZ];
    double m_divAccelField[GRID_NX][GRID_NY][GRID_NZ];
    double m_accelFieldCurr[3][GRID_NX][GRID_NY][GRID_NZ];
    double m_accelFieldFiltered[3][GRID_NX][GRID_NY][GRID_NZ];
    double m_pressureField[GRID_NX][GRID_NY][GRID_NZ];
    double m_RHS[GRID_NX][GRID_NY][GRID_NZ];

    size_t m_splinesCount;
    size_t m_timeStepCount;
    vector<ParticleTrack> m_vel;   // per-track velocity samples; only populated by loaders this build no longer ships (see README "Notes")
    vector<ParticleTrack> m_accel; // per-track acceleration samples; same caveat as m_vel
    std::vector<string> m_fileNames;
};

#endif // FILELOADER_H