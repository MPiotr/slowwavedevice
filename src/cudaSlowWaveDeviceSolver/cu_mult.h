#ifndef __CU_MULT_H
#define __CU_MULT_H
//#include "cuda_runtime.h"
//#include <shrUtils.h>
//#include <shrQATest.h>
//#include <cutil_inline.h>
extern cudaPitchedPtr  d2_rJ3, d2_iJ3, d2_int_rJ3, d2_int_iJ3, d2_W;
extern cudaPitchedPtr  d1_rJ3, d1_iJ3, d1_int_rJ3, d1_int_iJ3, d1_W;

/* int Nz = 0;
 double Lsolver = 0;
 double MultiplierGroupSpeedCoefficient = 0;*/

extern __device__  void biAverage(double *A, double *B, int p0, int datasize, int logsize);
#endif 
