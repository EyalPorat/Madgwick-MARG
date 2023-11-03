/*
MIT License

Copyright (c) 2023 Eyal Porat

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#pragma once

#include "quaternions.h"
#include <stdint.h>
#include <string.h>
#include <math.h>


struct MARG
{
	MARG(float beta_, float zeta_);
	
	//the filter step
	void updateMARG(double Dt, float mag[3], float gyro[3], float acceleration[3], float total_angle[3]);
	
  float total_ang[3]; //holding the calculated total angle in degrees
	Quaternion SEq; //the final estimated quaternion of orientation in relation to the global frame
	
	float beta; //the filter gain
	float zeta; //gyroscope bias drift compensation gain
	
	
	
	private:
		Quaternion m; //the magnetic readings
		Quaternion w; //the angular rate readings
		Quaternion a; //the acceleration readings
	
		//the gradient descent variables
		float f[6];
		float J[6][4];
		
		
		//estimated gyroscope bias drift
		Quaternion w_error;
		Quaternion w_b;
		
		Quaternion b; //holding the estimated direction of flux
		Quaternion h; //the magnetic readings rotated to the global frame
		Quaternion SEqDot_omega; //the estimated angular rate quaternion
		Quaternion SEqDotHat; //the gradient descent output
		
		//pre-calculated to reduce calculations
		Quaternion halfSEq;
		Quaternion twoSEq;
		Quaternion twob;
		Quaternion twom;
		
		
	
};

void vectorRot(Quaternion effector, Quaternion input_vector, float output_vector[3]);
