/*
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef KFPSimdPaddedAllocator_H
#define KFPSimdPaddedAllocator_H
#include <limits>
#include <new>
#include <type_traits>


namespace KFP
{
  namespace SIMD
  {

    /// @brief Allocator that pads allocations to a multiple of the SIMD width.
    /// @tparam T
    /// @tparam Width
    /// @tparam Alignment
    template<class T, std::size_t Width, std::size_t Alignment>
    class PaddedAllocator {

      static_assert(Width > 0);
      static_assert(Alignment > 0 && (Alignment & (Alignment - 1)) == 0);

      static constexpr std::size_t actual_alignment = Alignment < alignof(T) ? alignof(T) : Alignment;

     public:
      using value_type      = T;
      using is_always_equal = std::true_type;

      template<class U>
      struct rebind {
        using other = PaddedAllocator<U, Width, Alignment>;
      };

      PaddedAllocator() noexcept = default;

      template<class U>
      PaddedAllocator(const PaddedAllocator<U, Width, Alignment>&) noexcept
      {
      }

      [[nodiscard]] T* allocate(std::size_t n)
      {
        constexpr auto max_count = std::numeric_limits<std::size_t>::max() / sizeof(T);

        const auto extra = (Width - n % Width) % Width;

        if (n > max_count || extra > max_count - n) { throw std::bad_array_new_length{}; }

        const auto padded_n = n + extra;

        return static_cast<T*>(::operator new(padded_n * sizeof(T), std::align_val_t{actual_alignment}));
      }

      void deallocate(T* p, std::size_t) noexcept { ::operator delete(p, std::align_val_t{actual_alignment}); }
    };

    template<class T, class U, std::size_t W, std::size_t A>
    bool operator==(const PaddedAllocator<T, W, A>&, const PaddedAllocator<U, W, A>&) noexcept
    { return true; }

    template<class T, class U, std::size_t W, std::size_t A>
    bool operator!=(const PaddedAllocator<T, W, A>&, const PaddedAllocator<U, W, A>&) noexcept
    { return false; }


  }  // namespace SIMD
}  // namespace KFP

#endif  // KFPSimdPaddedAllocator
