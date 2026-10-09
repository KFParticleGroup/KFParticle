// -*- C++ Header -*-
/*
==================================================
Authors: A.Mithran;
Emails: mithran@fias.uni-frankfurt.de
==================================================
*/

#ifndef SIMD_SSE_CONSTANTS_H
#define SIMD_SSE_CONSTANTS_H

#include <cstddef>
#include <experimental/simd>

namespace stdx = std::experimental;

namespace KFP
{
  namespace SIMD
  {

    constexpr int SimdSize{16};
    constexpr int SimdLen{4};

  }  // namespace SIMD
}  // namespace KFP

#endif  // !SIMD_SSE_CONSTANTS_H
