#pragma once

#include <cmath>
#include <vector>

#include "../ParticleStructures/MCParticleStructure.h"

// This function is just an example, but not
// a real analysis template. Please, do not
// use it for ALICE data analysis - it is
// just senseless.

// Since covariance matrix is a simmetric matrix, we
// can consider only low triangle of this matrix.
// Also we will decompose it to rows, and write
// this rows one by one into vector.
//
// x, y, z, px, py, pz - our global parameters
//
//  Cov. matrix:
//
// σxx
// σxy   σyy
// σxz   σyz    σzz
// σxpx  σypx   σzpx    σpxpx
// σxpy  σypy   σzpy    σpxpy   σpypy
// σxpz  σypz   σzpz    σpxpz   σpypz   σpzpz
//
// xx, yy, zz,... - in cm^2
// pxpx, pypy, ... - in GeV^2
// xpx, ypy, zpz, ... - in cm*GeV

std::vector<float> MakeCovMatrix(const MCParticleStruct &part) {
  // 6x6 matrix has 21 indep. elements. Lets define all of them.
  std::vector<float> covmat(21);

  float pt = -1.;
  if (part.isMother)
    pt = part.finalPT();
  else
    pt = part.initialPT();

  covmat[0] =
      -3.59281e-06 + 1.02589e-05 / (2.90359e-02 + exp(-1.23947 / pt)); // σxx
  covmat[1] = -1.03361e-04 - 2.13155e-05 / pt +
              8.35099e-05 * log(3.34901 + (1. / pt)); // σxy
  covmat[2] =
      -6.17021e-06 + 1.21346e-05 / (2.73856e-02 + exp(-1.23947 / pt)); // σyy
  covmat[3] = 4.75842e-06 - 2.16231e-07 / pt;                          // σxz
  covmat[4] = 5.55892e-06 - 3.34455e-07 / pt;                          // σyz
  covmat[5] =
      -1.24273e-04 + 2.09914e-04 / (1.94591e-01 + exp(-1.25084 / pt)); // σzz
  covmat[6] = -9.77974e-06 - 1.26268e-07 * pow(1. / pt, 3) -
              4.13286e-06 * log(1. / pt); // σxpx
  covmat[7] = 1.04643e-06 + 5.47427e-09 * pow(1. / pt, 3) +
              1.02030e-07 * log(1. / pt); // σypx
  covmat[8] = -6.07308e-07 - 6.06728e-11 * pow(1. / pt, 5) -
              2.71488e-07 * log(1. / pt); // σzpx
  covmat[9] =
      2.46840e-05 + 1.73881e-06 / (4.59998e-03 + exp(-2.28556 * pt)); // σpxpx
  covmat[10] = 1.44735e-06 + 5.22400e-09 * pow(1. / pt, 3) -
               1.47596e-09 * log(1. / pt); // σxpy
  covmat[11] = 2.43044e-04 - 9.99721e-08 * pow(1. / pt, 3) -
               8.62114e-05 * log(1.79634e+01 + 1. / pt); // σypy
  covmat[12] = -2.09525e-06 - 1.63962e-08 * pow(1. / pt, 3) -
               1.37700e-06 * log(1. / pt); // σzpy
  covmat[13] = 7.32129e-05 + 2.35547e-05 * pt -
               7.30680e-05 * log(2.63928e+00 + pt); // σpxpy
  covmat[14] =
      2.15703e-05 + 1.3404e-06 / (4.31348e-03 + exp(-2.29888 * pt)); // σpypy
  covmat[15] = 3.98892e-07 - 1.70505e-09 * pow(1. / pt, 3) +
               8.76380e-08 * log(1. / pt); // σxpz
  covmat[16] = 4.56639e-07 + 3.78791e-11 * pow(1. / pt, 3) +
               4.69908e-08 * log(1. / pt); // σypz
  covmat[17] = 4.60901e-04 - 1.04157e-07 * pow(1. / pt, 3) -
               1.85200e-04 * log(1.35021e+01 + 1. / pt); // σzpz
  covmat[18] =
      -5.45974e-07 - 3.80233e-07 * pow(pt, 3) + 2.88266e-07 * log(pt); // σpxpz
  covmat[19] = 3.19847e-04 - 3.22758e-04 * pow(1. / pt, 0.1) +
               3.55611e-05 * log(-4.93008e-02 + 1. / pt); // σpypz
  covmat[20] =
      1.57546e-05 + 2.39963e-06 / (9.59981e-03 + exp(-2.02381 * pt)); // σpzpz

  // go from cm to microns
  /* covmat[0] *= 1e8;
   covmat[1] *= 1e8;
   covmat[2] *= 1e8;
   covmat[3] *= 1e8;
   covmat[4] *= 1e8;
   covmat[5] *= 1e8;

   covmat[6] *= 1e4;
   covmat[7] *= 1e4;
   covmat[8] *= 1e4;
   covmat[10] *= 1e4;
   covmat[11] *= 1e4;
   covmat[12] *= 1e4;
   covmat[15] *= 1e4;
   covmat[16] *= 1e4;
   covmat[17] *= 1e4;*/

  return covmat;
}