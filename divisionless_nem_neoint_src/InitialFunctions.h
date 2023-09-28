#pragma once

#include <cmath>
#include <utility>

#include <FixVec.h>

namespace AMDiS
{
  inline auto rotatedEllipse(WorldVector<double> c0, double eps, double dim, double angle, double a = 1.0, double b = 1.0, bool dirichlet = false)
  {
    WorldVector<double> R;
    R[0] = std::cos(angle * M_PI/180.0);
    R[1] = std::sin(angle * M_PI/180.0);

    WorldVector<double> c;
    c[0] = R[0] * c0[0] - R[1] * c0[1];
    c[1] = R[1] * c0[0] + R[0] * c0[1];

    return [eps,dim,R_=std::move(R),c_=std::move(c),angle,a,b,dirichlet](WorldVector<double> const& x)
    {
      WorldVector<double> x_;

      x_[0] = R_[0] * x[0] - R_[1] * x[1];
      x_[1] = R_[1] * x[0] + R_[0] * x[1];

      double dist = 0.5*(a + b)*(std::sqrt( sqr((x_[0] - c_[0])/a) + sqr((x_[1] - c_[1])/b) ) - 1.0);
      double result = -std::tanh(dist / (std::sqrt(2.0) * eps));
      
      if (!dirichlet){
      // add up all periodic realizations
      for (int i = 0; i < 2; ++i) {
        for (int d : {-1,1}) {
          WorldVector<double> c0(c_); c0[i] += d*dim;

          double dist = 0.5*(a + b)*(std::sqrt( sqr((x_[0] - c0[0])/a) + sqr((x_[1] - c0[1])/b) ) - 1.0);
          result = std::max(result, -std::tanh(dist / (std::sqrt(2.0) * eps)));
        }
      }        

      }

      return result;
    };

  }

  inline auto invertedRotatedEllipse(WorldVector<double> c0, double eps, double dim, double angle, double a = 1.0, double b = 1.0, bool dirichlet = false)
  {
    WorldVector<double> R;
    R[0] = std::cos(angle * M_PI/180.0);
    R[1] = std::sin(angle * M_PI/180.0);

    WorldVector<double> c;
    c[0] = R[0] * c0[0] - R[1] * c0[1];
    c[1] = R[1] * c0[0] + R[0] * c0[1];

    return [eps,dim,R_=std::move(R),c_=std::move(c),angle,a,b,dirichlet](WorldVector<double> const& x)
    {
      WorldVector<double> x_;

      x_[0] = R_[0] * x[0] - R_[1] * x[1];
      x_[1] = R_[1] * x[0] + R_[0] * x[1];

      double dist = 0.5*(a + b)*(std::sqrt( sqr((x_[0] - c_[0])/a) + sqr((x_[1] - c_[1])/b) ) - 1.0);
      double result = std::tanh(dist / (std::sqrt(2.0) * eps));
      /*
      if (!dirichlet){
      // add up all periodic realizations
      for (int i = 0; i < 2; ++i) {
        for (int d : {-1,1}) {
          WorldVector<double> c0(c_); c0[i] += d*dim;

          double dist = 0.5*(a + b)*(std::sqrt( sqr((x_[0] - c0[0])/a) + sqr((x_[1] - c0[1])/b) ) - 1.0);
          result = std::max(result, std::tanh(dist / (std::sqrt(2.0) * eps)));
        }
      }        

      }*/

      return result;
    };

  }

  inline auto annularRing(WorldVector<double> c0, double eps, double dim, double angle, double a_in = 1.0, double gap = 0.1, bool dirichlet = false)
  {
    //a_in is radius of inner cell, gap is the distance to outer cell from inner cell. radius of outer cell = a_in + gap
    WorldVector<double> R;
    R[0] = std::cos(angle * M_PI/180.0);
    R[1] = std::sin(angle * M_PI/180.0);

    WorldVector<double> c;
    c[0] = R[0] * c0[0] - R[1] * c0[1];
    c[1] = R[1] * c0[0] + R[0] * c0[1];

    return [eps,dim,R_=std::move(R),c_=std::move(c),angle,a_in,gap,dirichlet](WorldVector<double> const& x)
    {
      WorldVector<double> x_;

      x_[0] = R_[0] * x[0] - R_[1] * x[1];
      x_[1] = R_[1] * x[0] + R_[0] * x[1];

      double dist = 0.0;
      double result = 0.0;
      if(std::sqrt( sqr(x_[0] - c_[0]) + sqr(x_[1] - c_[1]) ) < a_in + gap/2){
        dist = (std::sqrt( sqr(x_[0] - c_[0]) + sqr(x_[1] - c_[1]) ) - a_in);
        result = -std::tanh(dist / (std::sqrt(2.0) * eps));
      } else {
        dist = (std::sqrt( sqr(x_[0] - c_[0]) + sqr(x_[1] - c_[1]) ) - a_in - gap);
        result = std::tanh(dist / (std::sqrt(2.0) * eps));
      }
      return result;
    };

  }

  inline auto noCell() //to assign empty phase
  {
    double result;
    return[](WorldVector<double> const& x)
    {
      double result = -1.0;
      return result; 
    };
  }	

  inline auto rectangle(WorldVector<double> c0, double eps, double dim, double a = 1.0, double b = 1.0, bool dirichlet = false)
  {
    return [c0,eps,dim,a,b,dirichlet](WorldVector<double> const& x)
    {
      double dist = (std::max(std::abs(x[0] - c0[0]) - 0.5*a,std::abs(x[1] - c0[1])-0.5*b) );
      double result = -std::tanh(dist / (std::sqrt(2.0) * eps));

      if (!dirichlet){
      // add up all periodic realizations
      for (int i = 0; i < 2; ++i) {
        for (int d : {-1,1}) {
          WorldVector<double> c_(c0); c_[i] += d*dim;

          double dist = (std::max(std::abs(x[0] - c_[0]) - 0.5*a,std::abs(x[1] - c_[1])-0.5*b) );
          result = std::max(result, -std::tanh(dist / (std::sqrt(2.0) * eps)));
        }
      }
      }


      return result;
    };
  } 
    inline auto invertedRectangle(WorldVector<double> c0, double eps, double dim, double a = 1.0, double b = 1.0, bool dirichlet = false)
  {
    return [c0,eps,dim,a,b,dirichlet](WorldVector<double> const& x)
    {
      double dist = (std::max(std::abs(x[0] - c0[0]) - 0.5*a,std::abs(x[1] - c0[1])-0.5*b) );
      double result = std::tanh(dist / (std::sqrt(2.0) * eps));
      /*
      if (!dirichlet){
      // add up all periodic realizations
      for (int i = 0; i < 2; ++i) {
        for (int d : {-1,1}) {
          WorldVector<double> c_(c0); c_[i] += d*dim;

          double dist = (std::max(std::abs(x[0] - c_[0]) - 0.5*a,std::abs(x[1] - c_[1])-0.5*b) );
          result = std::max(result, std::tanh(dist / (std::sqrt(2.0) * eps)));
        }
      }
      }*/


      return result;
    };
  } 

  inline auto hexagon(WorldVector<double> c0, double eps, double dim, double r = 1.0)
  {
    return [c0,eps,dim,r](WorldVector<double> const& x)
    {
      // Helpers
      double k0 = -0.866025404;
      double k1 = 0.37;
      double k2 = 0.577350269;

      WorldVector<double> q;
      q[0] = std::abs(x[0]-c0[0]);
      q[1] = std::abs(x[1]-c0[1]);

      double kDotQ = k0 * q[1] + k1 * q[0];
      q[0] -= 2*std::min(kDotQ,0.0) * k1;
      q[1] -= 2*std::min(kDotQ,0.0) * k0;

      q[1] -= q[1] < (-1.0) * k2 ? (-1.0) * k2 : (q[1] > k2 ? k2 : q[1]);
      q[0] -= r;

      double dist = 0.0;
      if (q[0] < 0.0)
        dist = (-1.0) * std::sqrt(sqr(q[0]) + sqr(q[1]));
      else 
        dist = std::sqrt(sqr(q[0]) + sqr(q[1]));

      double result = -std::tanh(dist / (std::sqrt(2.0) * eps));

      // add up all periodic realizations
      for (int i = 0; i < 2; ++i) {
        for (int d : {-1,1}) {
          WorldVector<double> c_(c0); c_[i] += d*dim;

          WorldVector<double> q;
          q[0] = std::abs(x[0]-c_[0]);
          q[1] = std::abs(x[1]-c_[1]);

          double kDotQ = k0 * q[1] + k1 * q[0];
          q[0] -= 2*std::min(kDotQ,0.0) * k1;
          q[1] -= 2*std::min(kDotQ,0.0) * k0;

          q[1] -= q[1] < (-1.0) * k2 ? (-1.0) * k2 : (q[1] > k2 ? k2 : q[1]);
          q[0] -= r;

          double dist = 0.0;
          if (q[0] < 0.0)
            dist = (-1.0) * std::sqrt(sqr(q[0]) + sqr(q[1]));
          else 
            dist = std::sqrt(sqr(q[0]) + sqr(q[1]));

          result = std::max(result, -std::tanh(dist / (std::sqrt(2.0) * eps)));
        }
      }
      


      return result;
    };
  } 


  inline auto trapezoidal(WorldVector<double> c0, double eps, double dim, double a = 1.0, double b = 1.0, double fac = 1.0, double angle = 0.0, bool dirichlet = false)
  {
    WorldVector<double> R;
    R[0] = std::cos(angle * M_PI/180.0);
    R[1] = std::sin(angle * M_PI/180.0);

    WorldVector<double> c;
    c[0] = R[0] * c0[0] - R[1] * c0[1];
    c[1] = R[1] * c0[0] + R[0] * c0[1];

    return [c0_=std::move(c),eps,dim,a,b,fac,R_=std::move(R),dirichlet](WorldVector<double> const& x)
    {
      WorldVector<double> x_;

      x_[0] = R_[0] * x[0] - R_[1] * x[1];
      x_[1] = R_[1] * x[0] + R_[0] * x[1];

      double dist = (std::max(std::abs(x_[0] - c0_[0]) - fac * (x_[1] - c0_[1]) - 0.5*a,std::abs(x_[1] - c0_[1])-0.5*b) );
      double result = -std::tanh(dist / (std::sqrt(2.0) * eps));


      return result;
    };
  } 

  inline auto arc(WorldVector<double> c0, double eps, double ta = 0.0, double tb = 0.0, double ra = 0.0, double rb = 0.0)
  {
    WorldVector<double> sca,scb;
    sca[0] = std::sin(ta); // Rotation counter-clockwise
    sca[1] = std::cos(ta);

    scb[0] = std::sin(tb); // arclength
    scb[1] = std::cos(tb);

    return [c0,eps,sca,scb,ra,rb](WorldVector<double> const& x)
    {
      WorldVector<double> x_; //rescale to new center
      x_[0] = x[0] - c0[0];
      x_[1] = x[1] - c0[1];

      WorldVector<double> p;
      p[0] = std::abs(sca[0] * x_[0] - sca[1] * x_[1]);
      p[1] = sca[1] * x_[0] + sca[0] * x_[1];

      double k = 0.0;
      if (scb[1] * p[0] <= scb[0] * p[1])
        k = std::sqrt(sqr(p[0]) + sqr(p[1]));

      double dist = std::sqrt(sqr(p[0]) + sqr(p[1]) + ra * ra - 2.0 * ra * k) - rb;
      double result = -std::tanh(dist / (std::sqrt(2.0) * eps));

      return result;
    };
  }   

} // end namespace AMDiS
