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

#include "MARG.h"

#define sq(x) ((x)*(x))

//calculate constants
const float pi = 3.141592654f;
const float rad_to_deg = 180 / pi;



MARG::MARG(float beta_, float zeta_)
{
	memset(f, 0, 6*sizeof(float));
	memset(J, 0, 4*6*sizeof(float));
	
	beta = beta_;
	zeta = zeta_;
	
	SEq.c[0] = 1;

}


//rotate a vector using a quaternion
void vectorRot(Quaternion effector, Quaternion input_vector, float output_vector[3])
{
	Quaternion output;
	output = effector * input_vector * effector.conj();
	output_vector[0] = output.c[1];
	output_vector[1] = output.c[2];
	output_vector[2] = output.c[3];
}


void MARG::updateMARG(double Dt, float mag[3], float gyro[3], float acceleration[3], float total_angle[3])
{
	
	//fill in the accel data and normalize it
	a.c[1] = acceleration[0];
	a.c[2] = acceleration[1];
	a.c[3] = acceleration[2];
	
	a.normalize();
	
	//fill in the gyro data and convert it to rad/s
	w.c[1] = gyro[0] / rad_to_deg;
	w.c[2] = gyro[1] / rad_to_deg;
	w.c[3] = gyro[2] / rad_to_deg;
	
	//fill in the mag data and normalize is
	m.c[1] = mag[0];
	m.c[2] = mag[1];
	m.c[3] = mag[2];
	
	m.normalize();
	
	
	//batch compute future used variables
  twoSEq.c[0] = 2.0f * SEq.c[0];
  twoSEq.c[1] = 2.0f * SEq.c[1];
  twoSEq.c[2] = 2.0f * SEq.c[2];
  twoSEq.c[3] = 2.0f * SEq.c[3];

  twob.c[1] = 2.0f * b.c[1];
  twob.c[2] = 2.0f * b.c[2];
  twob.c[3] = 2.0f * b.c[3];
	
	//compute the gradient descent step
	
	f[0] = 2 * (SEq.c[1] * SEq.c[3] - SEq.c[0] * SEq.c[2]) - a.c[1];
	f[1] = 2 * (SEq.c[0] * SEq.c[1] + SEq.c[2] * SEq.c[3]) - a.c[2];
	f[2] = 2 * (0.5f - sq(SEq.c[1]) - sq(SEq.c[2])) - a.c[3];
	
	f[3] = twob.c[1] * (0.5f - sq(SEq.c[2]) - sq(SEq.c[3])) + twob.c[3] * (SEq.c[1] * SEq.c[3] - SEq.c[0] * SEq.c[2]) - m.c[1];
	f[4] = twob.c[1] * (SEq.c[1] * SEq.c[2] - SEq.c[0] * SEq.c[3]) + twob.c[3] * (SEq.c[0] * SEq.c[1] + SEq.c[2] * SEq.c[3]) - m.c[2];
	f[5] = twob.c[1] * (SEq.c[0] * SEq.c[2] + SEq.c[1] * SEq.c[3]) + twob.c[3] * (0.5f - sq(SEq.c[1]) - sq(SEq.c[2])) - m.c[3];
	
	
	J[0][0] = -twoSEq.c[2];
	J[0][1] = twoSEq.c[3];
	J[0][2] = -twoSEq.c[0];
	J[0][3] = twoSEq.c[1];
	J[1][0] = twoSEq.c[1];
	J[1][1] = twoSEq.c[0];
	J[1][2] = twoSEq.c[3];
	J[1][3] = twoSEq.c[2];
	J[2][0] = 0;
	J[2][1] = -4.0f * SEq.c[1];
	J[2][2] = -4.0f * SEq.c[2];
	J[2][3] = 0;
	
	J[3][0] = -twob.c[3] * SEq.c[2];
	J[3][1] = twob.c[3] * SEq.c[3];
	J[3][2] = -4 * b.c[1] * SEq.c[2] - twob.c[3] * SEq.c[0];
	J[3][3] = -4 * b.c[1] * SEq.c[3] + twob.c[3] * SEq.c[1];
	J[4][0] = -twob.c[1] * SEq.c[3] + twob.c[3] * SEq.c[1];
	J[4][1] = twob.c[1] * SEq.c[2] + twob.c[3] * SEq.c[0];
	J[4][2] = twob.c[1] * SEq.c[1] + twob.c[3] * SEq.c[3];
	J[4][3] = -twob.c[1] * SEq.c[0] + twob.c[3] * SEq.c[2];
	J[5][0] = twob.c[1] * SEq.c[2];
	J[5][1] = twob.c[1] * SEq.c[3] - 4 * b.c[3] * SEq.c[1];
	J[5][2] = twob.c[1] * SEq.c[0] - 4 * b.c[3] * SEq.c[2];
	J[5][3] = twob.c[1] * SEq.c[1];
	
	//compute Jt * f
	
	SEqDotHat.c[0] = J[0][0] * f[0] + J[1][0] * f[1] + J[2][0] * f[2] + J[3][0] * f[3] + J[4][0] * f[4] + J[5][0] * f[5] + 1e-7;
	SEqDotHat.c[1] = J[0][1] * f[0] + J[1][1] * f[1] + J[2][1] * f[2] + J[3][1] * f[3] + J[4][1] * f[4] + J[5][1] * f[5] + 1e-7;
	SEqDotHat.c[2] = J[0][2] * f[0] + J[1][2] * f[1] + J[2][2] * f[2] + J[3][2] * f[3] + J[4][2] * f[4] + J[5][2] * f[5] + 1e-7;
	SEqDotHat.c[3] = J[0][3] * f[0] + J[1][3] * f[1] + J[2][3] * f[2] + J[3][3] * f[3] + J[4][3] * f[4] + J[5][3] * f[5] + 1e-7;
	
	//normalize the quaternion
	
	SEqDotHat.normalize();
	
	//calculate the gyroscope error direction
	
	w_error = twoSEq.conj() * SEqDotHat;
	
	
	w_b.c[1] += w_error.c[1] * Dt * zeta;
	w_b.c[2] += w_error.c[2] * Dt * zeta;
	w_b.c[3] += w_error.c[3] * Dt * zeta;
	
	w.c[1] -= w_b.c[1];
	w.c[2] -= w_b.c[2];
	w.c[3] -= w_b.c[3];
	
	//rotate the previous estimation by the new w to get SEq dot omega
	
	SEqDot_omega = (SEq * 0.5f) * w;
	
	//fuse the different parts
	
	SEqDotHat = SEqDotHat * beta;
	
	SEqDot_omega = SEqDot_omega - SEqDotHat;
	
	SEq = SEq + SEqDot_omega * Dt;
	
	SEq.normalize();
	
	//translate quaternion to euler angles
	
	total_ang[0] = atan2f(2 * (SEq.c[2] * SEq.c[3] + SEq.c[0] * SEq.c[1]), 1 - 2 * (sq(SEq.c[1]) + sq(SEq.c[2])));
	
	total_ang[1] = 2.0f * (SEq.c[0] * SEq.c[2] - SEq.c[3] * SEq.c[1]);
	total_ang[1] = (total_ang[1] >= 1.0f ? 1.0f : total_ang[1]);
	total_ang[1] = (total_ang[1] <= -1.0f ? -1.0f : total_ang[1]);
	total_ang[1] = asinf(total_ang[1]);	
	
	total_ang[2] = atan2f(2 * (SEq.c[1] * SEq.c[2] + SEq.c[0] * SEq.c[3]), 1 - 2 * (sq(SEq.c[2]) + sq(SEq.c[3])));

  total_ang[0] *= rad_to_deg;
  total_ang[1] *= rad_to_deg;
  total_ang[2] *= rad_to_deg;
	
	total_angle[0] = total_ang[0];
	total_angle[1] = total_ang[1];
	total_angle[2] = total_ang[2];
	
	
	//calculate the estimated flux direction
	h = SEq * m * SEq.conj();

  b.c[1] = sqrt(sq(h.c[1]) + sq(h.c[2]));
  b.c[3] = h.c[3];
	
}
