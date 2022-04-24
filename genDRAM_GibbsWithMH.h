//
//  genDRAM_GibbsWithMH.h
//
//
//  Created by X.Z on 10/4/18.
//

#ifndef genDRAM_GibbsWithMH_h
#define genDRAM_GibbsWithMH_h

#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
//#include <Eigen/core>
#include <random>
#include <numeric>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <cstdlib>
#include <fstream>
#include <string>
#include <chrono>
#include "PurineSynthesis.h"
#include "data_process.h"
//#include "Constants.h"
//#include "some_extern.h"


extern vector<int> ifFit;
extern vector<double> sd_tune;
extern    double mmYobs [3][10], mmYobs2 [3][10];
extern    double Yode [10],Yode2[10] ;
extern    double startp [124];
extern    double muprior [62], muprior2 [62];
extern    double sigmaprior [62];
//extern double beta_tmp[95], alpha_tmp[95];
extern    int observed [10];
extern double beta_tmp[20], alpha_tmp[20];
extern double beta_err[20], alpha_err[20];

extern const double t0;
extern const double tf;
extern const double df;

extern const int nruns;
extern  const  int nburnin;
extern const int nthinning;

extern const int N_tmp;
extern const int n_Y;
extern const int n_err;
extern const int n_err2;

void lpostLik(double& lpost, double *ptmp, double *err_curr, int myindex, int diffterm);

void methodGibbs(int& jruns);

template <typename Derived>
void reloCov(int index, Eigen::MatrixBase<Derived>& CAdapt);

void genDRAM_GibbsWithMH(int randseed, int diffterm);


#endif /* genDRAM_GibbsWithMH_h */
