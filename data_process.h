/*
 * data_process.h
 *
 *  Created on: Nov 25, 2019
 *      Author: zhangxu
 */

#ifndef DATA_PROCESS_H_
#define DATA_PROCESS_H_


#include <stdio.h>
#include <boost/array.hpp>

using namespace std;

typedef boost::array< double , 512 > state_type;

double* data_process(double data_arr[],state_type &Yd);


#endif /* DATA_PROCESS_H_ */
