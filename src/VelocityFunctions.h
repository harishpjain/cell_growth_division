#pragma once

#include <cmath>
#include <utility>

#include <FixVec.h>

namespace AMDiS
{
  inline auto elongateVelocity(WorldVector<double> center_of_mass_, double v0_, double major_axis_angle_, double scale_)
  {
    return [center_of_mass_, v0_, major_axis_angle_, scale_](WorldVector<double> const& x)
    {
      double x_elon_ = ((x[0]-center_of_mass_[0]) * std::cos(major_axis_angle_)) + ((x[1]-center_of_mass_[1])*std::sin(major_axis_angle_));
      //x_elon_[1] = R_[1] * x[0] + R_[0] * x[1]; not needed

      double result = (v0_ + x_elon_*scale_);
      return result;
    };

  }
}