#pragma once

namespace AMDiS
{

  /// \brief A Range implementation based on `first` and `last` iterators. The increment and
  /// dereference operations can be modified using functors.
  template <class Iter, class Range>
  struct IndexMap
  {
    using value_type = typename Range::value_type;

    struct const_iterator
        : public std::iterator<std::forward_iterator_tag, IndexMap::value_type>
    {
      using value_type = IndexMap::value_type;

      const_iterator(Iter it, Range const& range)
        : it_(it)
        , range_(range)
      {}

      const_iterator& operator++()
      {
        ++it_;
        return *this;
      }
      const_iterator  operator++(int) { const_iterator tmp(*this); ++(*this); return tmp; }

      bool operator==(const_iterator const& other) const { return it_ == other.it_; }
      bool operator!=(const_iterator const& other) const { return !(*this == other); }

      value_type operator*() const { return range_[*it_]; }

      Iter it_;  // the current iterator
      Range const& range_; // a range of values
    };

    /// \brief Constructor, stores the begin and end iterator `first` and `last`, respectively
    /// and a range with values.
    /**
      * \param first An iterator pointing at the begin of an index sequence
      * \param last  An iterator pointing past the end of an index sequence
      * \param range A range that is accessed indirectly by indices of the index sequence
      **/
    IndexMap(Iter first, Iter last, Range const& range)
      : first_(first)
      , last_(last)
      , range_(range)
    {}

    /// Provide an iterator to the begin of the range
    const_iterator cbegin() const { return {first_, range_}; }
    /// Provide an iterator to the end of the range
    const_iterator cend()   const { return {last_, range_}; }

    /// Provide an iterator to the begin of the range
    const_iterator begin() const { return cbegin(); }
    /// Provide an iterator to the end of the range
    const_iterator end()   const { return cend(); }

    /// Returns the size of the range, i.e. the distance of begin and end iterator
    std::size_t size() const { return std::distance(first_, last_); }


  private:

    Iter first_;
    Iter last_;
    Range range_;
  };


  /// \brief Generator for \ref IndexMap, \relates IndexMap
  /**
    * Requirement:
    * - `Iter` models `ForwardIterator`
    *
    **/
  template <class Iter, class Range>
  IndexMap<Iter, Range>
  index_map(Iter first, Iter last, Range const& range)
  {
    return {first, last, range};
  }

} // end namespace AMDiS
