#include "fileloader.h"
#include <filesystem>

ParticleTrackLoader::ParticleTrack::ParticleTrack()
{
}

ParticleTrackLoader::ParticleTrack::ParticleTrack(double ix, double iy, double iz, int it)
{
    x.push_back(ix);
    y.push_back(iy);
    z.push_back(iz);
    t.push_back(it);
}

ParticleTrackLoader::ParticleTrackLoader(string iPath)
{
    m_path = iPath;
    for (size_t i = 0; i < 10000000; i++)
    {
        m_numbers[i] = -1;
    }
}

// Lists every "*.dat" file directly inside m_path, sorted alphabetically (this
// is what turns "particles0.dat, particles1.dat, ..." into an ordered sequence
// of time steps). Originally used a POSIX dirent.h + qsort/VLA combination that
// doesn't build with MSVC; replaced with std::filesystem + std::sort, same result.
// message) instead of crashing -- the original called readdir() on a null
// DIR* returned by a failed opendir(), which is undefined behavior.
void ParticleTrackLoader::findDataFiles()
{
    namespace fs = std::filesystem;
    m_fileNames.clear();

    std::error_code dirError;
    fs::directory_iterator dirIt(m_path, dirError);
    if (dirError)
    {
        printf("Could not open data directory \"%s\" (%s) -- starting with no tracks loaded\n",
               m_path.c_str(), dirError.message().c_str());
        return;
    }

    for (const auto& entry : dirIt)
    {
        if (entry.path().extension() == ".dat")
        {
            m_fileNames.push_back(entry.path().string());
        }
    }
    std::sort(m_fileNames.begin(), m_fileNames.end());
    for (size_t i = 0; i < m_fileNames.size(); i++)
    {
        printf("FilesName= %s\n", m_fileNames[i].c_str());
    }
}

// Reads every file in m_fileNames (one per time step, in order) and appends each
// line's (x,y,z) sample to the track identified by the trailing particle-id
// column, building m_data (one entry per distinct particle id seen).
void ParticleTrackLoader::loadParticleTracks()
{
    m_data.clear();
    m_timeStepCount = m_fileNames.size();
    for (size_t i = 0; i < m_fileNames.size(); i++)
    {
        FILE *file_data = fopen(m_fileNames[i].c_str(), "r");
        printf("name=%s\n", m_fileNames[i].c_str());
        if (file_data == nullptr)
        {
            printf("Could not open track file \"%s\" -- skipping it\n", m_fileNames[i].c_str());
            continue;
        }
        char str[128];
        if (fgets(str, sizeof(str), file_data) == nullptr ||
            fgets(str, sizeof(str), file_data) == nullptr ||
            fgets(str, sizeof(str), file_data) == nullptr)
        {
            fclose(file_data);
            continue;
        }
        double xx, yy, zz;
        int num;
        while (fscanf(file_data, "%lf %lf %lf %d", &xx, &yy, &zz, &num) == 4)
        {
            if ((num >= 0) && (num < 1000000))
            {
                if (m_numbers[num] == -1)
                {
                    m_numbers[num] = (int)m_data.size();
                    m_data.push_back(ParticleTrack(xx, yy, zz, (int)i));
                }
                else
                {
                    m_data[m_numbers[num]].x.push_back(xx);
                    m_data[m_numbers[num]].y.push_back(yy);
                    m_data[m_numbers[num]].z.push_back(zz);
                    m_data[m_numbers[num]].t.push_back((int)i);
                    xmin = xx < xmin ? xx : xmin;
                    xmax = xx > xmax ? xx : xmax;
                    ymin = yy < ymin ? yy : ymin;
                    ymax = yy > ymax ? yy : ymax;
                    zmin = zz < zmin ? zz : zmin;
                    zmax = zz > zmax ? zz : zmax;
                }
            }
        }
        fclose(file_data);
    }
    m_splinesCount = m_data.size();
}