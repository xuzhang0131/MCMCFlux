/*
 * PurineSynthesis.h
 *
 *  Created on: Oct 15, 2018
 *      Author: zhangxu
 */

#ifndef PURINESYNTHESIS_H_
#define PURINESYNTHESIS_H_

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include <time.h>
//#include "some_extern.h"
//#include "Constants.h"

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double , 512 > state_type;

extern const double MoleculeNumberInOneNanoMole;// = 6.02214129e14;
extern const int NumberOfIsotopomers;

extern double odeParameters[57];

//typedef boost::array< double , 608 > state_type;

extern const int CO2;
extern const int NH3;
extern const int cCO2;
extern const int cNH3;
extern const int mCO2;
extern const int mNH3;
extern const double AvogadroConstant;

//extern double odeParameters[57];

void PurineSynthesis( const state_type &Xd , state_type &dxdt , double t);


#endif /* PURINESYNTHESIS_H_ */
