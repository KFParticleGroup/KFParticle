// -*- C++ Header -*-
/*
==================================================
Authors: A.Mithran;
Emails: mithran@fias.uni-frankfurt.de
==================================================
*/

#ifndef SIMD_SSE_TYPE_H
#define SIMD_SSE_TYPE_H

#include "float32.h"
#include "int32.h"
#include "mask32.h"

#include <stdexcept>
#include <type_traits>

namespace KFP
{
  namespace SIMD
  {
    using float_v = Float32_128;
    static_assert(std::is_same<float_v::value_type, float>::value,
                  "[Error]: Invalid value type for SSE float SimdClass.");

    using int_v = Int32_128;
    static_assert(std::is_same<int_v::value_type, int>::value, "[Error]: Invalid value type for SSE int SimdClass.");

    using float_m = Mask32_128;
    using int_m   = Mask32_128;

    KFP_SIMD_INLINE Int32_128 toInt(const Float32_128& a) { return _mm_cvtps_epi32(a.simd()); }
    KFP_SIMD_INLINE Float32_128 toFloat(const Int32_128& a) { return _mm_cvtepi32_ps(a.simd()); }

    KFP_SIMD_INLINE Int32_128 reinterpretAsInt(const Float32_128& a) { return _mm_castps_si128(a.simd()); }
    KFP_SIMD_INLINE Float32_128 reinterpretAsFloat(const Int32_128& a) { return _mm_castsi128_ps(a.simd()); }

    inline const int_v gkIndicesSequenceI(int_v::indicesSequenceTmp());

    inline const float_v gkIndicesSequenceF(toFloat(gkIndicesSequenceI));

    KFP_SIMD_INLINE float_v gather(const float_v::value_type* data, const int_v& indices)
    {
      float_v result;
      result.gatherTmp(data, indices);
      return result;
    }

    KFP_SIMD_INLINE int_v gather(const int_v::value_type* data, const int_v& indices)
    {
      int_v result;
      result.gatherTmp(data, indices);
      return result;
    }

    template<int N>
    KFP_SIMD_INLINE float_v rotate(const float_v& v)
    { return v.rotateTmp<N>(); }

    template<int N>
    KFP_SIMD_INLINE int_v rotate(const int_v& v)
    { return v.rotateTmp<N>(); }

  }  // namespace SIMD
}  // namespace KFP

#endif  // !SIMD_SSE_TYPE_H
