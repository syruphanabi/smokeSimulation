#ifndef _MYWORLD_
#define _MYWORLD_

#include <vector>

#define IX(i, j) ((i)+(getNumCells()+2)*(j))
#define SWAPPING(x0,x) {double *tmp=x0;x0=x;x=tmp;}

class MyWorld {
 public:
    bool Cshape;
    bool Fan;

    MyWorld();

    virtual ~MyWorld();

    void initialize(int _numCells, double _timeStep, double _diffCoef, double _viscCoef);
    double getTimeStep();
    int getNumCells() { return mNumCells;}
    double getDensityR(int _index) { return mDensityR[_index]; }
    double getDensityG(int _index) { return mDensityG[_index]; }
    double getDensityB(int _index) { return mDensityB[_index]; }
    double getVelocityU(int _index) { return mU[_index]; }
    double getVelocityV(int _index) { return mV[_index]; }
    void setDensityR(int _i, int _j, double _source) { mDensityR[IX(_i, _j)] += mTimeStep * _source; }
    void setDensityG(int _i, int _j, double _source) { mDensityG[IX(_i, _j)] += mTimeStep * _source; }
    void setDensityB(int _i, int _j, double _source) { mDensityB[IX(_i, _j)] += mTimeStep * _source; }
    void setU(int _i, int _j, double _force) { mU[IX(_i, _j)] += mTimeStep * _force; }
    void setV(int _i, int _j, double _force) { mV[IX(_i, _j)] += mTimeStep * _force; }

    void simulate();

 protected:
    void densityStep(double *_x, double *_x0);
    void velocityStep(double *_u, double *_v, double *_u0, double *_v0);
    void diffuseDensity(double *_x, double *_x0);
    void diffuseVelocity(double *_u, double *_v, double *_u0, double *_v0);
    void advectDensity(double *_d, double *_d0, double *_u, double *_v);
    void advectVelocity(double *_u, double *_v, double *_u0, double *_v0);
    void project(double *_u, double *_v, double *_u0, double *_v0);
    void externalForces();
    void linearSolve(double *_x, double *_x0, double _a, double _c);
    void setBoundary(double *_x);
    void setVelocityBoundary(double *_u, double *_v);
    void setBlock(double *_x);

    int mNumCells;
    double mTimeStep;
    double mDiffusionCoef;
    double mViscosityCoef;
    double *mU;
    double *mV;
    double *mPreU;
    double *mPreV;
    // double *mDensity;
    double *mDensityR;
    double *mDensityG;
    double *mDensityB;
    // double *mPreDensity;
    double *mPreDensityR;
    double *mPreDensityG;
    double *mPreDensityB;
};

#endif
