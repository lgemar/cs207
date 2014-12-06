#ifndef SPACE_SEARCHER_HPP
#define SPACE_SEARCHER_HPP
/** @file Space.hpp
 * @brief Define the SpaceSearcher class for making spatial searches.
 */

#include <cassert>

#include "Point.hpp"
#include "BoundingBox.hpp"
#include "MortonCoder.hpp"

/** @class SpaceSearcher
 * @brief Class for making spatial searches that uses MortonCoder as backend.
 */
template <typename T, typename T2Point, int L = 5>
class SpaceSearcher {
public:

  ////////////////////////////////////
  // TYPE DEFINITIONS AND CONSTANTS //
  ////////////////////////////////////

  /** Type of indexes and sizes. */
  typedef unsigned size_type;

  /** The number of levels in the MortonCoder. */
  static constexpr int NumLevels = 5;
  /** Type of MortonCoder. */
  typedef MortonCoder<NumLevels> MortonCoderType;
  /** Type of the Morton codes. */
  typedef typename MortonCoderType::code_type code_type;

  /** Type of iterators, which iterate over items inside a BoundingBox. */
  class Iterator;
  /** Synonym for Iterator */
  typedef Iterator iterator;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Constructor. */
  template <typename TIter>
  SpaceSearcher(TIter t_begin, TIter t_end, T2Point t2p)
      : mc_(), c2t_(), t2p_(t2p) {

    // For performance checks
#if 0
    CS207::Clock clock;
    clock.start();
    std::cout << "\nConstructing SpaceSearcher..." << std::endl;
#endif
    // Determine pmin and pmax, i.e., the minimum and maximum corners
    // for the bounding box.
    Point pmin = t2p_(*t_begin);
    Point pmax = pmin;
    for (auto it = t_begin; it != t_end; ++it) {
      Point p = t2p_(*it);
      for (int i = 0; i < 3; ++i) {
        if      (p[i] < pmin[i]) pmin[i] = p[i];
        else if (p[i] > pmax[i]) pmax[i] = p[i];
      }
    }

    // Create MortonCoder instance.
	bb_ = BoundingBox(pmin, pmax);
    mc_ = MortonCoderType(bb_);

    // Map Morton codes to data items.
    c2t_.clear();
    c2t_.resize(mc_.end_code);
    for (auto it = t_begin; it != t_end; ++it) {
      Point p = t2p_(*it);
      c2t_[mc_.code(p)].push_back(*it);
    }

    // Sanity check
    size_type item_count = 0;
    size_type max_items = 0;
    for (auto it = c2t_.begin(); it != c2t_.end(); ++it) {
      auto& items = *it;
      item_count += items.size();
      if (items.size() > max_items)
        max_items = items.size();
    }
#if 0
    std::cout << "Time: " << clock.seconds() << " seconds." << std::endl;
    std::cout << "Total number of elements = " << item_count << std::endl;
    std::cout << "Total number of cells = " << mc_.end_code << std::endl;
    std::cout << "Max. number of elements per cell = " << max_items << std::endl;
    std::cout << std::endl;
#endif
  }

  /** Default destructor */
  ~SpaceSearcher() = default;

  //////////////
  // ITERATOR //
  //////////////

  /** @class SpaceSearcher::Iterator
   * @brief Iterator class for data items. A forward iterator. */
  class Iterator : private totally_ordered<Iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef T value_type;
    /** Type of pointers to elements. */
    typedef T* pointer;
    /** Type of references to elements. */
    typedef T& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid Iterator. */
    Iterator() : s_(nullptr), bb_(), code_(-1), loc_(0) {
    }
    /** Method to dereference an iterator.
     * @pre This is a valid Iterator.
     */
    T operator*() const {
      assert(is_valid());
      return s_->c2t_[code_][loc_];
    }
    /** Method to increment an iterator.
     *
     * Note that the return value may be end(), and therefore invalid.
     */
    Iterator& operator++() {
      assert(s_ != nullptr && !bb_.empty());
      ++loc_;
      while (code_ < s_->mc_.code(bb_.max())+1) {
        while (loc_ < s_->c2t_[code_].size()) {
          // Make sure that item really is inside of BoundingBox.
          if (bb_.contains(s_->t2p_(s_->c2t_[code_][loc_])))
            return *this;
          ++loc_;
        }
        // Advance to next Morton cell that overlaps with BoundingBox.
        code_ = s_->mc_.advance_to_box(code_+1, s_->mc_.code(bb_.min()),
                                                s_->mc_.code(bb_.max()));
        loc_ = 0;
      }

      // If no more valid Morton codes, we return end().
      code_ = s_->mc_.code(bb_.max())+1;
      loc_ = 0;
      return *this;
    }
    /** Method to compare two iterators. */
    bool operator==(const Iterator& x) const {
      return ((s_ == x.s_) &&
              (bb_.min() == x.bb_.min()) &&
              (bb_.max() == x.bb_.max()) &&
              (code_ == x.code_) &&
              (loc_ == x.loc_));
    }

   private:
    // Allow SpaceSearcher to access Iterator's private members.
    friend class SpaceSearcher;
    // Pointer back to the SpaceSearcher container.
    SpaceSearcher* s_;
    // BoundingBox associated with this Iterator.
    BoundingBox bb_;
    // Morton code of current item.
    code_type code_;
    // Index of current item inside its Morton cell
    // (note that different items may have the same Morton code).
    size_type loc_;
    /** Private constructor. */
    Iterator(const SpaceSearcher* s, BoundingBox bb, code_type code, size_type loc)
        : s_(const_cast<SpaceSearcher*>(s)), bb_(bb), code_(code), loc_(loc) {
      // Advance Iterator to valid position, if necessary.
      fix();
    }
    /** Helper method to determine if this Iterator is valid.
     *
     * Note that we don't actually check that the point is contained inside
     * the BoundingBox. This is handled by fix() and operator++().
     */
    bool is_valid() const {
      if (s_ == nullptr || bb_.empty())
        return false;
      if (code_ >= s_->mc_.code(bb_.max())+1)
        return false;
      return loc_ < s_->c2t_[code_].size();
    }
    /** Helper method to advance this Iterator until it reaches a valid
     * position, or end().
     */
    void fix() {
      assert(s_ != nullptr && !bb_.empty());
      if (code_ >= s_->mc_.code(bb_.max())+1) {
        // Make equal to end() and return.
        code_ = s_->mc_.code(bb_.max())+1;
        loc_ = 0;
        return;
      }
      if (loc_ >= s_->c2t_[code_].size() ||
          !bb_.contains(s_->t2p_(s_->c2t_[code_][loc_]))) {
        // Advance to next valid item (or end()).
        operator++();
      }
    }
  };

  /** Method to graph the bounding box of this space */
  BoundingBox bounding_box() {
  	return bb_;
  }

  /** Method to return an iterator pointing to the first item
   * in a given BoundingBox.
   */
  iterator begin(const BoundingBox& bb) const {
    assert(!bb.empty());
    return Iterator(this, bb, mc_.code(bb.min()), 0);
  }

  /** Method to return an iterator pointing to "one past"
   * the last item in a given BoundingBox.
   */
  iterator end(const BoundingBox& bb) const {
    assert(!bb.empty());
    return Iterator(this, bb, mc_.code(bb.max())+1, 0);
  }

private:
  // MortonCoder instance associated with this SpaceSearcher.
  MortonCoderType mc_;

  // Bounding box associated with this spatial searcher
  BoundingBox bb_;

  // Structure that maps Morton codes to data items of type T.
  std::vector<std::vector<T>> c2t_;

  // Mapping from data items to points.
  T2Point t2p_;
};

#endif
