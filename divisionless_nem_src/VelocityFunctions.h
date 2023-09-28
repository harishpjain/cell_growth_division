#pragma once

#include <cmath>
#include <utility>

#include <FixVec.h>

namespace AMDiS
{
  inline auto elongateVelocity(double xscale_, double yscale_, double dscale_)
  {
    return [xscale_, yscale_, dscale_](WorldVector<double> const& x, WorldVector<double> center_of_mass_, double v0_, double theta_, double radius)
    {
      //X_elon is the distance from center of mass of cell along the direction of motion. 
      //y_elon is the distance from line of motion of cell measured relative to center of mass 
      double x_elon_ = ((x[0]-center_of_mass_[0]) * std::cos(theta_)) + ((x[1]-center_of_mass_[1])*std::sin(theta_));
      double y_elon_ = ((x[0]-center_of_mass_[0]) * std::cos(theta_+(M_PI/2))) + ((x[1]-center_of_mass_[1])*std::sin(theta_+(M_PI/2)));
      double d_elon_ = 0.0;//((x[0]-center_of_mass_[0]) * std::cos(theta_+(3*M_PI/4))) + ((x[1]-center_of_mass_[1])*std::sin(theta_+(3*M_PI/4)));
      double elon_ = ((x[0]-center_of_mass_[0]) * std::cos(theta_+M_PI)) + ((x[1]-center_of_mass_[1])*std::sin(theta_+M_PI));
      //x_elon_[1] = R_[1] * x[0] + R_[0] * x[1]; not needed

      //double result = (v0_ - std::fabs(x_elon_)*scale_);
      if(elon_>1.5*radius) elon_ = 1.5*radius; //softmax function to elon_ to limit elongation.
      if(elon_<-1.5*radius) elon_ = -1.5*radius;
      //double result = v0_*(1 + x_elon_*xscale_ - std::fabs(y_elon_)*yscale_ + d_elon_*dscale_ + elon_*dscale_);
      double result = v0_*(1+elon_*dscale_);
      return result;
    };

  }

  inline auto elongateVelocityAlt(WorldVector<double> center_of_mass_, double v0_, double theta_, double scale_) 
  //alternative function that should do the same as above
  {
    return [center_of_mass_, v0_, theta_, scale_](WorldVector<double> const& x)
    {
      double b = std::atan2(x[1]-center_of_mass_[1], x[0]-center_of_mass_[0]);
      double a = theta_ - b;
      double c = std::cos(a);
      double norm_x = std::sqrt(std::pow(x[0]-center_of_mass_[0], 2) + std::pow(x[1]-center_of_mass_[1], 2));
      double x_elon_ = norm_x * c;
      //x_elon_[1] = R_[1] * x[0] + R_[0] * x[1]; not needed

      //double result = (v0_ - std::fabs(x_elon_)*scale_);
      double result = v0_ + x_elon_*scale_;
      return result;
    };

  }

  inline auto shearY(double domainDimensiony_)
  {
    return [domainDimensiony_](WorldVector<double> const& x, double v0_)
    {
      return v0_*std::abs((domainDimensiony_/2.0 - x[1]))/(domainDimensiony_/2.0);
    };

  }
}