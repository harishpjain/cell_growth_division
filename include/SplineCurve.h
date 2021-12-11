#pragma once

namespace AMDiS
{
  class SplineCurve
  {
  public:
    template <class Points>
    SplineCurve(Points const& points)
      : points_(points.begin(), points.end())
    {
      parametrization_.resize(points.size());
      adjacent_distance(points.begin(), points.end(), parametrization_.begin(),
                        [](WorldVector<double> const& p0, WorldVector<double> const& p1) { return norm(p1 - p0); });
      
      // normalize to [0,1]
      double length = parametrization_.back();
      std::transform(parametrization_.begin(), parametrization_.end(), parametrization_.begin(),
                     [length](double p) { return p/length; });
    }
    
    WorldVector<double> operator()(double t) const
    {
      auto it = std::lower_bound(parametrization_.begin(), parametrization_.end(), t);
      if (it == parametrization_.begin())
        return points_[0];
      
      std::size_t i = std::distance(parametrization_.begin(), it);
      
      double lambda = (t - parametrization_[i-1])/(parametrization_[i] - parametrization_[i-1]);
      return lambda * points_[i] + (1.0-lambda) * points_[i-1];
    }
        
  private:
    
    template <class InputIt, class OutputIt, class BinaryOperation>
    OutputIt adjacent_distance(InputIt first, InputIt last, 
                               OutputIt d_first, BinaryOperation op)
    {
        if (first == last) return d_first;
    
        typedef typename std::iterator_traits<InputIt>::value_type value_t;
        value_t acc = *first;
        *d_first = 0;
        while (++first != last) {
            value_t val = *first;
            *++d_first = op(val, acc);
            acc = std::move(val);
        }
        return ++d_first;
    }
        
  private:
    
    std::vector<WorldVector<double>> points_;
    std::vector<double> parametrization_;
  };

} // end namespace AMDiS
