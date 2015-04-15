/*
 * egsheader.h
 *
 *  Created on: Apr 15, 2015
 *      Author: boris
 */

#ifndef EGSHEADER_H_
#define EGSHEADER_H_

#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <complex>
#include <stdio.h>
#include <assert.h>
#include <cstdlib>
#include <assert.h>
#include <stdint.h>
#include <cmath>
//#include <tgmath.h>

#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include "sys/time.h"
#define _USE_MATH_DEFINES
#include <time.h>
#include<sys/stat.h>
#include<sys/types.h>
#include <vector>


#include <omp.h>

using namespace std;
//								Commands Definitions
#define createHostVariable(type, name,n) type *name; name  = new type[n]()
#define DeleteHostVariable(name) delete[] name

void sayhello();

template <typename T>
void MemoryCopy(T *in, T *out, int nx, int ny);

template <typename T>
void PorosityFunction( T *Phi, T *iPhi,T *Pw, T *Sxx, T *Szz, T phiini, T phir,  T phia, int nx, int ny);

template <typename T>
void PermeabilitySinglePhase(T *Kappa, T *Phi, T *Kappa0,T *Pw, T *Sxx, T *Szz,  T phiini, T cap_depth, T k_cap, T c,T h, T e, T kappa_hole,  T dx, T dy,  int nx, int ny);

template <typename T>
void ThermalConductivityFunction(T *TCm, T *Frac, T *Rcp, T *Phi, T *Rhow, T *Cpw, T tcr, T rhor0, T cpr, T  tcw, int nx, int ny);

template <typename T>
void StateEqTemperetureDensityTransport(T *Viscw , T *Rhof, T *Tw, T *Pw, T *Cs, T *Betaw, T *Cpw, T rhof0, int nx, int ny);

template <typename T>
void LTEDensityTransportSML(T *Tr, T *iTr, T *Tw,T *Vwx,  T *Vwy, T *Phi,  T *Rhow, T*Rcp, T *Cpw, T *TCr, T *x,T *y, T h, T e,T tw_fresh, T dx, T dy, T t, T dt,int nx, int ny);

template <typename T>
void LTNENonIsoThermal(T *Tw, T *Tr,  T *iTw, T *iTr, T *Vwx, T *Vwy, T *Phi, T *Rhow, T *Cpw, T tcw,  T rhor0, T  cpr, T tcr, T rockr, T *x, T *y,T tw_fresh, T dx, T dy, T dt, int nx, int ny);

template <typename T>
void LTNENonIsoThermalTR(T *iTw, T *Tr,  T *iTr, T *Phi, T *Hta, T  t, T rhor0, T  cpr,  T tcr,T h,T e,T tw_fresh, T dx, T dy, T dt, int nx, int ny);

template <typename T>
void LTNENonIsoThermalTF(T *Tw, T *iTw, T *iTr,T *Vwx,T *Vwy,T *Phi,T *Rhow,T *Hta,T *Cpw,T t,T tcw,T *x,T *y,T h,T e,T tw_fresh, T dx,T dy,T dt,int nx,int ny);

template <typename T>
void Hfs(T *Hta ,T *TCr,T *Phi, T *Re, T tcw,T rockr ,int nx, int ny);

template <typename T>
void EGSPressureFunction(T *Pw,T *Vwx, T *Vwy,T *iPw, T *Rhof,T *Phi, T *Kappa, T *Viscw ,T *Betaf, T *Betaphi, T *S, T *Re, T *Pw0, T *Y, T rockr , T h, T e,T pmax, T pmin, T vwx0 ,T dx, T dy, T t, int g_flag, T dt,int nx, int ny);

template <typename T>
void EGSMassFraction(T *Cs, T *iCs, T *Dm, T *Vwx, T *Vwy, T *Phi, T rhow0, T rho_max,T h, T e,T dx, T dy, T dt, int nx, int ny);

template <typename T>
void TimeStepCalculation(T &dt,T *Phi, T *Kappa, T *Viscw, T *Betaphi, T *Betaw, T *Rhow,T *TCm,T *Dm,T *Cpw,T tcr,T tcw,T cpr,T rhor0,T dx,T dy,int nx, int ny);

template <typename T>
 void AdvectionGamma( T *Gamma, T *Alpha, T *Rhof ,T *Cpw, T *Phi,T rhor0, T cpr,T dx, T dy, int nx, int ny);

 template <typename T>
void AdvectionCoefficient( T *Adv, T *Rhof ,T *Cpw, T *Phi,T rhor0, T cpr, int lte_flag,int nx, int ny);

template <typename T>
void Fracture(T *K, T *Phi,T *TCr, T *Pw, T *Szz, T *Sxx, T k_fracture,  T perm_max,T phi_fracture,T depth_fracture, T width_fracture,T theta, T fracture_aperture, T l_fracture, T tcr_fracture,T phiini, T phir, T phia, T c, double PI, T dx, T dy, int nx, int ny);

template <typename T>
void SaveDistAve(T *M0, T *M1, T *M2,T *M3,T *M4,T *M5, T *M6, T *M7, T *M8, string name, T dx,T dy,  int tim, int nx,int ny);

template <typename T>
void SaveVector(T *Mx,  T *My, string name, T dx, T dy, int tim, int nx, int ny);

template <typename T>
void SaveVtk(T *Pw,T *Cs,T *Tf, T *Tr, T *Kappa, T *Phi,T *Vx,  T *Vy, string name, T dx, T dy, int tim, int nx, int ny);

template <typename T>
void DensityTransportMrkrCll(T *TCs,T *iCs,T *Vwx,T *Vwy,T *Xm, T *Ym,T xmin,T ymin,T h,T e,T marknum, T dx,T dy,T dt,int nx,int ny);

template <typename T>
void DensityTransportSML(T *Cs,T *iCs,T *Dm,T *Phi,T *Vwx,T *Vwy,T *x,T *y, T * Rhof,T h,T e, T cini,T dx,T dy,T t, T teq,T dt,int nx,int ny);

template <typename T>
void DensityTransportWENO5RK(T *Cs,T *iCs,T *Dm,T *Vwx,T *Vwy,T h, T e,T dx,T dy,T dt,int nx,int ny);

void indexsearch5(int *index, int i, int indi, int indf);

template <typename T>
T weno5FluxM(T *iCs,T *V,int flag, T maxa);

template <typename T>
T weno5FluxP(T *iCs,T *V,int flag, T maxa);

template <typename T>
void InterpolationPoints(int *xind,int *yind,T &npox,T &npoy,T *x, T *y,int nx,int ny);

template <typename T>
T spline_cubic(const T a[4],T x);

template <typename T>
T set2limits( T value, T min, T max );

template <typename T>
void BottomAverageFunction(T &ave, T *M, int nx, int ny);

template <typename T>
void TopAverageFunction( T &ave,  T *M ,int nx, int ny);

template <typename T>
void LeftAverageFunction( T &ave,  T *M ,int nx, int ny);

template <typename T>
void RightAverageFunction( T &ave,  T *M ,int nx, int ny);

template <typename T>
void MidAverageFunction(T &ave, T *M, int nx, int ny);

template <typename T>
void AverageFunction(T &ave, T *M, int nx, int ny);

template <typename T>
void InjectionAverageFunction(T &ave, T *M, T h, T e, T dy, int nx, int ny);

template <typename T>
void InjectionAverageFunction(T &ave, T *M, T h, T e, T dy, int nx, int ny);


template <typename T>
void readdata(T *odata,int nx,int ny, const char* s);


template <typename T>
void writedata(T *idata,int nx,int ny, int tim, string name);

template <typename T>
void PdWrite(T *Pd, T *Pw, T *Pw0,T d ,int nx,int ny, T dy, int tim, string name);

template <typename T>
void PressureRelaxFunction(T *Pw,T *Vwx, T *Vwy,T *iPw, T *Rhof ,T *Phi, T *Kappa, T *Viscw ,T *Betaf, T *Betaphi, T *S, T vwx0 ,T dx, T dy, T t, int g_flag, T dt,int nx, int ny);


#endif /* EGSHEADER_H_ */
