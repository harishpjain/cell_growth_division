#pragma once

#include <cmath>
#include <utility>

#include <FixVec.h>
//inspired by Initial Functions.h, and from Nonomura
namespace AMDiS{
	inline auto halfCutOne(double eps, double angle, WorldVector<double> center)
	{
		WorldVector<double> omega;
		omega[0] = std::cos(angle*M_PI/180.0);
		omega[1] = std::sin(angle*M_PI/180.0);
		
		return [eps, angle, omega_=std::move(omega), center](WorldVector<double> const&x)	{
			WorldVector<double> x_diff_;
			
			x_diff_[0] = x[0] - center[0];
			x_diff_[1] = x[1] - center[1];
			
			//double result = 0.0;
			double dot_prod = x_diff_[0]*omega_[0] + x_diff_[1]*omega_[1];
			double result = 1 + std::tanh(dot_prod/eps);
			
			return result;
		};
	}
	inline auto halfCutTwo(double eps, double angle, WorldVector<double> center)
	{
		WorldVector<double> omega;
		omega[0] = std::cos(angle*M_PI/180.0);
		omega[1] = std::sin(angle*M_PI/180.0);
		
		return [eps, angle, omega_=std::move(omega), center](WorldVector<double> const&x)	{
			WorldVector<double> x_diff_;
			
			x_diff_[0] = x[0] - center[0];
			x_diff_[1] = x[1] - center[1];
			
			//double result = 0.0;
			double dot_prod = x_diff_[0]*omega_[0] + x_diff_[1]*omega_[1];
			double result = 1 - std::tanh(dot_prod/eps);
			
			return result;
		};
	}
	
	inline auto rectangleCut(double eps, double angle, WorldVector<double> center, double length){
		//2*eps: width of cut, angle is orientation of rectangle, length: length of rectangle
		WorldVector<double> omega;
		omega[0] = std::cos(angle*M_PI/180.0);
		omega[1] = std::sin(angle*M_PI/180.0);
 		
		return [eps, angle, omega_=std::move(omega), center, length](WorldVector<double> const&x)	{
			WorldVector<double> cut_rot;
			cut_rot[0] = (x[0]-center[0])*omega_[0]+(x[1]-center[1])*omega_[1];
			cut_rot[1] = -(x[0]-center[0])*omega_[1]+(x[1]-center[1])*omega_[0];
			if(std::max(std::abs(cut_rot[0])-0.5*12*eps, std::abs(cut_rot[1])-0.5*length) <= 0){
				return 0.0;
			}
			else if(std::max(std::abs(cut_rot[0])-0.5*12*eps, std::abs(cut_rot[1])-0.5*length) <= 0){
				return 0.5;
			}
			else{
			return 1.0;
			}
		};
	}
	
	
		inline auto tanhCut(double eps, double angle, WorldVector<double> center)
	{
		WorldVector<double> omega;
		omega[0] = std::cos(angle*M_PI/180.0);
		omega[1] = std::sin(angle*M_PI/180.0);
		
		return [eps, angle, omega_=std::move(omega), center](WorldVector<double> const&x)	{
			WorldVector<double> x_diff_;
			
			x_diff_[0] = x[0] - center[0];
			x_diff_[1] = x[1] - center[1];
			
			double dot_prod = x_diff_[0]*omega_[0] + x_diff_[1]*omega_[1];
			//needs to change
			double result = std::tanh(dot_prod/eps);
			
			return result;
		};
	}	
}//end namespace AMDiS