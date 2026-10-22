/*
==================================================
Authors: S.Gorbunov;
Emails: se.gorbunov@gsi.de
==================================================
*/

#ifndef KFPSIMD_H
#define KFPSIMD_H

#include "KFPSimdPaddedVector.h"

#include <experimental/simd>
#include <stdexcept>
#include <type_traits>

//namespace stdx = std::experimental;
namespace stdx
{
  using namespace std::experimental;
  using namespace std::experimental::__proposed;
}  // namespace stdx


namespace kfp
{
  namespace simd
  {

#define KFP_SIMD_INLINE inline __attribute__((always_inline))

    // TODO: try with SIMD size 16

    using float_v = stdx::fixed_size_simd<float, 4>;
    using int_v   = stdx::rebind_simd_t<int, float_v>;

    using float_m = float_v::mask_type;
    using int_m   = int_v::mask_type;

    constexpr int SimdLen{float_v::size()};   // number of elements in the SIMD float type
    constexpr int SimdSize{sizeof(float_v)};  // size in bytes of the SIMD float type

    static_assert(std::is_same<float_v::value_type, float>::value,
                  "[Error]: Invalid value type for SSE float SimdClass.");

    static_assert(std::is_same<int_v::value_type, int>::value, "[Error]: Invalid value type for SSE int SimdClass.");


    //--- Casting between different mask types

    template<class ToMask, class T, class Abi>
    [[nodiscard]] ToMask mask_cast(const stdx::simd_mask<T, Abi>& src) noexcept
    {
      using FromMask = stdx::simd_mask<T, Abi>;

      static_assert(stdx::is_simd_mask_v<ToMask>, "ToMask must be a simd_mask");
      static_assert(ToMask::size() == FromMask::size(), "Masks must have the same number of lanes");

      if constexpr (std::is_constructible_v<ToMask, const FromMask&>) { return ToMask(src); }
      else {
        constexpr auto src_align = stdx::memory_alignment_v<FromMask>;
        constexpr auto dst_align = stdx::memory_alignment_v<ToMask>;
        constexpr auto alignment = src_align > dst_align ? src_align : dst_align;

        alignas(alignment) bool lanes[FromMask::size()];
        src.copy_to(lanes, stdx::vector_aligned);
        return ToMask(lanes, stdx::vector_aligned);
      }
    }

    //--- precomputed index sequences for int and float vectors

    inline const int_v gkIndicesSequenceI([](auto lane) { return static_cast<int_v::value_type>(lane); });

    inline const float_v gkIndicesSequenceF([](auto lane) { return static_cast<float_v::value_type>(lane); });

    //--- select() function for std::experimental::simd

    template<typename T>
    struct type_identity {
      using type = T;
    };

    template<typename T>
    using type_identity_t = typename type_identity<T>::type;

    template<typename T, typename Abi>
    KFP_SIMD_INLINE stdx::simd<T, Abi> select(const stdx::simd_mask<T, Abi>& mask,
                                              const type_identity_t<stdx::simd<T, Abi>>& on_true,
                                              const type_identity_t<stdx::simd<T, Abi>>& on_false)
    {
      auto result               = on_false;
      stdx::where(mask, result) = on_true;
      return result;
    }

    //-- TODO: SG: replace with std::simd::unchecked_gather_from() when switching to c++26

    KFP_SIMD_INLINE float_v gather(const PaddedVector<float_v>& data, const int_v& indices)
    {
      return float_v([&](auto lane) -> float_v::value_type { return data[indices[lane]]; });
    }

    KFP_SIMD_INLINE int_v gather(const PaddedVector<int_v>& data, const int_v& indices)
    {
      return int_v([&](auto lane) -> int_v::value_type { return data[indices[lane]]; });
    }

    KFP_SIMD_INLINE float_v gather(const float* data, const int_v& indices)
    {
      return float_v([&](auto lane) -> float_v::value_type { return data[indices[lane]]; });
    }

    KFP_SIMD_INLINE int_v gather(const int* data, const int_v& indices)
    {
      return int_v([&](auto lane) -> int_v::value_type { return data[indices[lane]]; });
    }

    //--- rotate() function for SIMD vectors

    template<int N, typename T, typename Abi>
    KFP_SIMD_INLINE stdx::simd<T, Abi> rotate(const stdx::simd<T, Abi>& v)
    {
      using V = stdx::simd<T, Abi>;

      constexpr int len   = static_cast<int>(V::size());
      constexpr int rem   = N % len;
      constexpr int shift = rem < 0 ? rem + len : rem;

      if constexpr (shift == 0) { return v; }
      else {
        return V([&v](auto i) -> T {
          constexpr int dst = static_cast<int>(decltype(i)::value);
          constexpr int src = (dst + len - shift) % len;
          return v[src];
        });
      }
    }

  }  // namespace simd
}  // namespace kfp

#endif  // !KFPSIMD_H
