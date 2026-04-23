#pragma once

#include <cmath>
#include <vector>

#include "../ParticleStructures/MCParticleStructure.h"
#include "ALICETrackParamConverter.h"

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
  double glob_covmat[21];
  double xyz[3];
  double pxpypz[3];
  short sign = part.charge;
  double qOverPt;

  if (part.isMother) {
    xyz[0] = part.finalX;
    xyz[1] = part.finalY;
    xyz[2] = part.finalZ;
    qOverPt = part.charge / part.finalPT();
  } else {
    pxpypz[0] = part.finalPX;
    pxpypz[1] = part.finalPY;
    pxpypz[2] = part.finalPZ;
    qOverPt = part.charge / part.initialPT();
  }

  ALICETrackParamsConverter glob_loc_track(xyz, pxpypz, glob_covmat, sign);

  double loc_x;
  glob_loc_track.GetX(loc_x);
  double loc_alpha;
  glob_loc_track.GetAlpha(loc_alpha);
  double loc_params[5];
  glob_loc_track.GetLocalParams(loc_params);
  double loc_covmat[15];
  loc_covmat[0] = ....;
  loc_covmat[1] = 0;
  loc_covmat[2] = ....;
  loc_covmat[3] = ....;
  loc_covmat[4] = 0;
  loc_covmat[5] = ....;
  loc_covmat[6] = 0;
  loc_covmat[7] = ....;
  loc_covmat[8] = 0;
  loc_covmat[9] = ....;
  loc_covmat[10] = ....;
  loc_covmat[11] = 0;
  loc_covmat[12] = ....;
  loc_covmat[13] = 0;
  loc_covmat[14] = ....;

  ALICETrackParamsConverter loc_glob_track(loc_params, loc_covmat, loc_x,
                                           loc_alpha, sign);

  loc_glob_track.GetGlobalCovMat(glob_covmat);

  std::vector<float> covmat(glob_covmat, glob_covmat + 21);
  return covmat;
}