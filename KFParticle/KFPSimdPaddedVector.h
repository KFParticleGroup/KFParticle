/*
==================================================
Authors: S.Gorbunov;
Emails: se.gorbunov@gsi.de
==================================================
*/

#ifndef KFPSIMD_PADDED_VECTOR_H
#define KFPSIMD_PADDED_VECTOR_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <experimental/simd>
#include <vector>

namespace kfp
{
  namespace simd
  {

    /// @brief SIMD-aligned amd padded vector of scalars - supports efficient access and manipulation of elements in SIMD blocks
    /// @tparam Simd The SIMD type used for alignment and padding. Simd::value_type is the scalar type of the elements.
    template<class Simd>
    class PaddedVector {
     public:
      using value_type = typename Simd::value_type;

      explicit PaddedVector(std::size_t n = 0, value_type value = {}) : fBlocks(blocksForSize(n), Simd(value)), fSize(n)
      {
      }

      // Copy constructor
      PaddedVector(const PaddedVector&) = default;

      // Copy assignment operator
      PaddedVector& operator=(const PaddedVector&) = default;

      // Move constructor. Keeps the source consistent
      PaddedVector(PaddedVector&& other) noexcept
        : fBlocks(std::move(other.fBlocks))
        , fSize(std::exchange(other.fSize, 0))
        , fIsPaddingFilled(std::exchange(other.fIsPaddingFilled, false))
      { other.fBlocks.clear(); }

      // Move assignment operator. Keeps the source consistent
      PaddedVector& operator=(PaddedVector&& other) noexcept
      {
        if (this != &other) {
          fBlocks          = std::move(other.fBlocks);
          fSize            = std::exchange(other.fSize, 0);
          fIsPaddingFilled = std::exchange(other.fIsPaddingFilled, false);
          other.fBlocks.clear();
        }
        return *this;
      }

      /// @brief Get the number of elements in the vector.
      /// @return The number of elements in the vector.
      std::size_t size() const noexcept { return fSize; }

      /// @brief Get the number of elements in a SIMD block.
      /// @return The number of elements in a SIMD block.
      static constexpr std::size_t width() noexcept { return fkWidth; }


      /// @brief Check if the vector is empty.
      /// @return True if the vector is empty, false otherwise.
      bool empty() const noexcept { return fSize == 0; }


      /// @brief Return the SIMD element proxy without converting it to a scalar.
      /// @param i The index of the element.
      /// @return The SIMD element proxy at the specified index.
      decltype(auto) operator[](std::size_t i)
      {
        assert(i < fSize);
        return fBlocks[i / fkWidth][i % fkWidth];
      }

      /// @brief Access an element by index.
      /// @param i The index of the element.
      /// @return The scalar element at the specified index.
      value_type operator[](std::size_t i) const
      {
        assert(i < fSize);
        return fBlocks[i / fkWidth][i % fkWidth];
      }

      /// @brief Get the number of SIMD blocks in the vector.
      /// @return The number of SIMD blocks in the vector.
      std::size_t blockCount() const noexcept { return fBlocks.size(); }

      /// @brief Access a SIMD block by index.
      /// @param i The index of the SIMD block.
      /// @return The SIMD block at the specified index.
      Simd& block(std::size_t i)
      {
        assert(i < fBlocks.size());
        return fBlocks[i];
      }

      /// @brief Access a SIMD block by index (const version).
      /// @param i The index of the SIMD block.
      /// @return The SIMD block at the specified index (const reference).
      const Simd& block(std::size_t i) const
      {
        assert(i < fBlocks.size());
        return fBlocks[i];
      }

      /// @brief Preserve existing elements and initialize new elements with value.
      /// @param n The new size of the vector.
      /// @param value The value to initialize new elements with (default is value-initialized).
      /// @details All lanes of newly allocated blocks are initialized.
      ///          Padding is not automatically updated to repeat the last element.
      void resize(std::size_t n, value_type value = {})
      {
        const auto old_size = fSize;

        fBlocks.resize(blocksForSize(n), Simd(value));

        // New elements within the previous partial block need explicit
        // initialization: vector::resize only initializes new blocks.
        if (n > old_size && old_size % fkWidth != 0) {
          const auto first_lane = old_size % fkWidth;
          const auto count      = std::min(n - old_size, fkWidth - first_lane);

          auto& previous_last = fBlocks[old_size / fkWidth];

          for (std::size_t lane = first_lane; lane < first_lane + count; ++lane) { previous_last[lane] = value; }
        }

        fSize            = n;
        fIsPaddingFilled = false;
      }

      /// @brief Fill unused lanes of the last block with the last logical element.
      /// @details Do nothing if the array is empty or the last block is full.
      void fillPaddingWithLastElement()
      {
        const auto used = fSize % fkWidth;

        if (used == 0) {
          fIsPaddingFilled = true;
          return;
        }

        auto& last             = fBlocks.back();
        const value_type value = last[used - 1];

        for (std::size_t lane = used; lane < fkWidth; ++lane) { last[lane] = value; }
        fIsPaddingFilled = true;
      }

      bool isPaddingFilled() const { return fIsPaddingFilled; }

      void markPaddingFilled() { fIsPaddingFilled = true; }

      /// @brief Compute the number of SIMD blocks needed for n elements.
      /// @param n The number of elements.
      /// @return The number of SIMD blocks required to store n elements.
      /// @details This computes ceil(n / fkWidth) without overflowing n + fkWidth - 1.
      static constexpr std::size_t blocksForSize(std::size_t n) noexcept { return n / fkWidth + (n % fkWidth != 0); }

     private:
      static constexpr std::size_t fkWidth{Simd::size()};  ///< The number of elements in a SIMD block.

      std::vector<Simd> fBlocks{};  ///< The SIMD blocks storing the elements.
      std::size_t fSize{0};         ///< The number of logical elements in the vector.

      bool fIsPaddingFilled{
        false};  ///< Indicates whether the padding of the last block has been filled with reasonable values (typically the last logical element).
    };

  }  // namespace simd
}  // namespace kfp


#endif  // !KFPSIMD_PADDED_VECTOR_H