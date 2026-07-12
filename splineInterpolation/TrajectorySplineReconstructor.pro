TEMPLATE = app
TARGET = TrajectorySplineReconstructor
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -lopenGL32 -lGLU32
SOURCES += \
        fileloader.cpp \
        globals.cpp \
        leastsquaressolver.cpp \
        main.cpp \
        splineinterpolator.cpp \
        tools.cpp
LIBS += my_lib/glut32.lib

HEADERS += \
    fileloader.h \
    globals.h \
    leastsquaressolver.h \
    splineinterpolator.h \
    tools.h
