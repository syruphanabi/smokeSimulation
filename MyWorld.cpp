#include "MyWorld.h"
#include <iostream>
#include <math.h>

using namespace std;


MyWorld::MyWorld() {
}

MyWorld::~MyWorld() {
}

void MyWorld::initialize(int _numCells, double _timeStep, double _diffCoef, double _viscCoef) {
    mNumCells = _numCells;
    mTimeStep = _timeStep;
    mDiffusionCoef = _diffCoef;
    mViscosityCoef = _viscCoef;

    int size = (mNumCells + 2) * (mNumCells + 2);

    // Allocate memory for velocity and density fields
    mU = new double[size];
    mV = new double[size];
    mPreU = new double[size];
    mPreV = new double[size];
    // mDensity = new double[size];
    // mPreDensity = new double[size];
    mDensityR = new double[size];
    mPreDensityR = new double[size];
    mDensityG = new double[size];
    mPreDensityG = new double[size];
    mDensityB = new double[size];
    mPreDensityB = new double[size];

    for (int i = 0; i < size; i++) {
        mU[i] = mPreU[i] = 0.0;
        mV[i] = mPreV[i] = 0.0;
        // mDensity[i] = mPreDensity[i] = 0.0;
        mDensityR[i] = mPreDensityR[i] = 0.0;
        mDensityG[i] = mPreDensityG[i] = 0.0;
        mDensityB[i] = mPreDensityB[i] = 0.0;
    }

    Cshape = false;
    Fan = true;
}

double MyWorld::getTimeStep() {
    return mTimeStep;
}

void MyWorld::simulate() {
    velocityStep(mU, mV, mPreU, mPreV);
    //cout << 'v' << endl;
    densityStep(mDensityR, mPreDensityR);
    densityStep(mDensityG, mPreDensityG);
    densityStep(mDensityB, mPreDensityB);
    //cout << 'd' << endl;
    if (Fan) externalForces();
    //cout << 'f' << endl;
}

void MyWorld::densityStep(double *_x, double *_x0) {
    SWAPPING(_x, _x0); // _x now points at mPreDensity
    diffuseDensity(_x, _x0); // Diffusion on _x which pointst at mPreDensity
    SWAPPING(_x, _x0); // _x now points at mDensity
    advectDensity(_x, _x0, mU, mV); // Advection on _x which points at mDensity
}

void MyWorld::velocityStep(double *_u, double *_v, double *_u0, double *_v0) {
    SWAPPING(_u, _u0); // _u now points at mPreU
    SWAPPING(_v, _v0); // _v now points at mPreV
    diffuseVelocity(_u, _v, _u0, _v0);
    project (_u, _v, _u0, _v0);
    SWAPPING(_u, _u0); // _u now points at mU
    SWAPPING(_v, _v0); // _u now points at mV
    advectVelocity(_u, _v, _u0, _v0);
    project(_u, _v, _u0, _v0);
}

void MyWorld::diffuseDensity(double *_x, double *_x0) {
    double a = mTimeStep * mDiffusionCoef * mNumCells * mNumCells;
    linearSolve(_x, _x0, a, 1 + 4 * a);
    setBoundary(_x);
}

void MyWorld::diffuseVelocity(double *_u, double *_v, double *_u0, double *_v0) {
    double a = mTimeStep * mViscosityCoef * mNumCells * mNumCells;
    linearSolve(_u, _u0, a, 1 + 4 * a);
    linearSolve(_v, _v0, a, 1 + 4 * a);
    setVelocityBoundary(_u, _v);
}

void MyWorld::advectDensity(double *_d, double *_d0, double *_u, double *_v) {
    double dt0 = mTimeStep * mNumCells;  // h / x
    for (int i = 1; i <= mNumCells; i++) {
        for (int j = 1; j <= mNumCells; j++) {
          //          if (abs(_u[IX(i, j)]) > 0.001)

            double x = i- dt0 * _u[IX(i,j)];  // dt0 * _u[IX(i,j)] computes how many cells can a particle travel in one time step
            double y = j - dt0 * _v[IX(i,j)];
            if (x < 0.5)
                x = 0.5f;
            if (x > mNumCells + 0.5)
                x = mNumCells + 0.5;
            int i0 = (int)x;
            int i1 = i0 + 1;
            if (y < 0.5)
                y = 0.5;
            if (y > mNumCells + 0.5)
                y = mNumCells + 0.5;
            int j0 = (int)y;
            int j1 = j0 + 1;
            double s1 = x - i0;
            double s0 = 1 - s1;
            double t1 = y - j0;
            double t0 = 1 - t1;
            _d[IX(i,j)] = s0 * (t0 * _d0[IX(i0, j0)] + t1 * _d0[IX(i0,j1)])+ s1 * (t0 * _d0[IX(i1, j0)] + t1 * _d0[IX(i1,j1)]);
	     }
    }
    setBoundary(_d);
}

void MyWorld::advectVelocity(double *_u, double *_v, double *_u0, double *_v0) {
    // TODO: Add velocity advection code here
    advectDensity(_u, _u0, _u0, _v0);
    advectDensity(_v, _v0, _u0, _v0);
}

void MyWorld::project(double *_u, double *_v, double *_u0, double *_v0) {
   // TODO: Add projection code here

   int i, j, k;
   double h = 1.0/mNumCells;
   int size = (mNumCells + 2) * (mNumCells + 2);

   double* mDiv = new double[size];
   double* mP = new double[size];

   for (i = 1; i <= mNumCells; i++){
     for (j = 1; j <= mNumCells; j++){
       mDiv[IX(i,j)] = - 0.5 * h * (_u0[IX(i+1,j)] - _u0[IX(i-1,j)] + _v0[IX(i, j+1)] - _v0[IX(i, j-1)]);
       mP[IX(i,j)] = 0.0;
     }
   }
   setBoundary(mDiv); setBoundary(mP);

   for (k = 0; k < 20; k++){
     for (i = 1; i <= mNumCells; i++){
       for (j = 1; j <= mNumCells; j++){
         mP[IX(i,j)] = (mDiv[IX(i,j)] + mP[IX(i+1,j)] + mP[IX(i-1,j)] + mP[IX(i,j+1)] + mP[IX(i,j-1)])/4;
       }
     }
     setBoundary(mP);
   }

   for (i = 1; i <= mNumCells; i++){
     for (j = 1; j <= mNumCells; j++){
       _u[IX(i,j)] -= 0.5 * (mP[IX(i+1,j)] - mP[IX(i-1,j)])/h;
       _v[IX(i,j)] -= 0.5 * (mP[IX(i,j+1)] - mP[IX(i,j-1)])/h;
     }
   }
   setVelocityBoundary(_u, _v);

}

void MyWorld::externalForces() {
    int size = (mNumCells + 2) * (mNumCells + 2);
    double norm;
    for (int i = 0; i< size; i++) {
        int x = i%(mNumCells+2);
        int y = (y == 0 ? i/(mNumCells+2)-1 : i/(mNumCells+2));
        double h = 0.0005;
        norm = sqrt(max(pow((x - mNumCells/2),2) + pow((y - mNumCells/2),2),1.0));
        if (norm > 5){
            if (Cshape){
                if (norm <= 15 or x >= 45){
                  mU[i] += (x - mNumCells/2)/norm * h - (y - mNumCells/2)/norm * 0.8 * h;
                  mV[i] += (y - mNumCells/2)/norm * h + (x - mNumCells/2)/norm * 0.8 * h;
              }
            }else{
              mU[i] += (x - mNumCells/2)/norm * h - (y - mNumCells/2)/norm * 0.8 * h;
              mV[i] += (y - mNumCells/2)/norm * h + (x - mNumCells/2)/norm * 0.8 * h;
            }
        }
    }
    //setVelocityBoundary(mU, mV);
}

void MyWorld::linearSolve(double *_x, double *_x0, double _a, double _c) {
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= mNumCells; i++) {
            for (int j = 1; j <= mNumCells; j++) {
                _x[IX(i, j)] = (_x0[IX(i, j)] + _a * (_x[IX(i-1, j)] + _x[IX(i+1, j)] + _x[IX(i, j-1)] + _x[IX(i, j+1)])) / _c;
            }
        }
    }
}

void MyWorld::setBoundary(double *_x) {
    for (int i = 1; i <= mNumCells; i++) {
        _x[IX(0 ,i)] = _x[IX(1,i)];
        _x[IX(mNumCells+1, i)] = _x[IX(mNumCells, i)];
        _x[IX(i, 0)] = _x[IX(i, 1)];
        _x[IX(i, mNumCells+1)] = _x[IX(i, mNumCells)];

    }
    _x[IX(0, 0)] = 0.5 * (_x[IX(1, 0)] + _x[IX(0, 1)]);
    _x[IX(0, mNumCells+1)] = 0.5 * (_x[IX(1, mNumCells+1)] + _x[IX(0, mNumCells)]);
    _x[IX(mNumCells+1, 0)] = 0.5 * (_x[IX(mNumCells, 0)] + _x[IX(mNumCells+1, 1)]);
    _x[IX(mNumCells+1, mNumCells+1)] = 0.5 * (_x[IX(mNumCells, mNumCells+1)] + _x[IX(mNumCells+1, mNumCells)]);

    if (Cshape) setBlock(_x);
}

void MyWorld::setVelocityBoundary(double *_u, double *_v) {
    for (int i = 1; i <= mNumCells; i++) {
        _u[IX(0 ,i)] = -_u[IX(1,i)];
        _u[IX(mNumCells+1, i)] = -_u[IX(mNumCells, i)];
        _u[IX(i, 0)] = _u[IX(i, 1)];
        _u[IX(i, mNumCells+1)] = _u[IX(i, mNumCells)];

        _v[IX(0 ,i)] = _v[IX(1,i)];
        _v[IX(mNumCells+1, i)] = _v[IX(mNumCells, i)];
        _v[IX(i, 0)] = -_v[IX(i, 1)];
        _v[IX(i, mNumCells+1)] = -_v[IX(i, mNumCells)];
    }
    _u[IX(0, 0)] = 0.5 * (_u[IX(1, 0)] + _u[IX(0, 1)]);
    _u[IX(0, mNumCells+1)] = 0.5 * (_u[IX(1, mNumCells+1)] + _u[IX(0, mNumCells)]);
    _u[IX(mNumCells+1, 0)] = 0.5 * (_u[IX(mNumCells, 0)] + _u[IX(mNumCells+1, 1)]);
    _u[IX(mNumCells+1, mNumCells+1)] = 0.5 * (_u[IX(mNumCells, mNumCells+1)] + _u[IX(mNumCells+1, mNumCells)]);
    _v[IX(0, 0)] = 0.5 * (_v[IX(1, 0)] + _v[IX(0, 1)]);
    _v[IX(0, mNumCells+1)] = 0.5 * (_v[IX(1, mNumCells+1)] + _v[IX(0, mNumCells)]);
    _v[IX(mNumCells+1, 0)] = 0.5 * (_v[IX(mNumCells, 0)] + _v[IX(mNumCells+1, 1)]);
    _v[IX(mNumCells+1, mNumCells+1)] = 0.5 * (_v[IX(mNumCells, mNumCells+1)] + _v[IX(mNumCells+1, mNumCells)]);

    if (Cshape) {setBlock(_u);setBlock(_v);}
}

void MyWorld::setBlock(double *_x){
    int size = (mNumCells + 2) * (mNumCells + 2);
    double norm;
    for (int i = 0; i< size; i++) {
        int x = i%(mNumCells+2);
        int y = (y == 0 ? i/(mNumCells+2)-1 : i/(mNumCells+2));
        norm = sqrt(max(pow((x - mNumCells/2),2) + pow((y - mNumCells/2),2),1.0));
        if (norm > 15 && norm < 20 && x < 45){
          _x[i] = 0;
        }
    }
}
