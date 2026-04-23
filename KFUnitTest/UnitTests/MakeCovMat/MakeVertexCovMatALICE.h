#pragma once
#define MAKE_VERTEX_COV_MAT

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
// x, y, z - our global parameters
//
//  Cov. matrix:
//
// σxx
// σxy   σyy
// σxz   σyz    σzz
//
// xx, yy, zz,... - in cm^2

std::vector<float> MakeVertexCovMatrix(const MCParticleStruct &part) {
  // 3x3 matrix has 6 indep. elements. Lets define all of them.
  std::vector<float> covmat(6);

  /*float pt = -1.;
  if (part.isMother)
      pt = part.finalPT();
  else
      pt = part.initialPT();*/

  covmat[0] = 1e-6;                                        // σxx
  covmat[1] = 1e-6 * (-1 + (double)rand() / RAND_MAX * 2); // σxy
  covmat[2] = 1e-6;                                        // σyy
  covmat[3] = 1e-6 * (-1 + (double)rand() / RAND_MAX * 2); // σxz
  covmat[4] = 1e-6 * (-1 + (double)rand() / RAND_MAX * 2); // σyz
  covmat[5] = 1e-6;                                        // σzz

  // Random number generator
  // std::random_device rd;
  // std::mt19937 gen(rd());
  // std::normal_distribution<double> dist(0.0, 1.0);

  // // Generate random values for covariance matrix elements
  // for (int i = 0; i < 6; ++i) {
  //     covmat[i] = dist(gen);
  // }

  // // Adjust elements to ensure positive semi-definiteness
  // covmat[0] = std::fabs(covmat[0]); // Ensure positive
  // covmat[2] = std::fabs(covmat[2]); // Ensure positive
  // covmat[5] = std::fabs(covmat[5]); // Ensure positive

  // // Ensure that off-diagonal elements are smaller than the square root of
  // the product of the corresponding diagonal elements covmat[1] =
  // std::min(std::fabs(covmat[1]), std::sqrt(covmat[0] * covmat[2])); covmat[3]
  // = std::min(std::fabs(covmat[3]), std::sqrt(covmat[0] * covmat[5]));
  // covmat[4] = std::min(std::fabs(covmat[4]), std::sqrt(covmat[2] *
  // covmat[5]));

  return covmat;
}