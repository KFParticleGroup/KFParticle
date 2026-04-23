#pragma once

#include <vector>

#include "../ParticleStructures/MCParticleStructure.h"
#include "ALICETrackParamConverter.h"

std::vector<float> MakeCovMatrix(const MCParticleStruct &part) {

  MCParticle(x, y, z, px, py, pz, ...)-- > -- > Convert from global to local-- >
      -- > {
    parameters
  }
  -- > -- > {
    PythonScript(parameters)-- > -- > RegressionTree(parameters)-- > -- >
        return CovMat[15] from PythonScript()
  }
  -- > -- > Convert from local to global-- > -- >
      return std::vector<float> covMat;
}