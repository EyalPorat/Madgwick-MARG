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

#include <math.h>
#include "quaternions.h"

#define sq(x) ((x)*(x))

Quaternion::Quaternion()
{
  c[0] = c[1] = c[2] = c[3] = 0;
}

Quaternion::Quaternion(float a0, float a1, float a2, float a3)
{
  c[0] = a0;
  c[1] = a1;
  c[2] = a2;
  c[3] = a3;
}

Quaternion Quaternion::conj() const
{
  return Quaternion(c[0], -1.00f * c[1], -1.00f * c[2], -1.00f * c[3]);
}

void Quaternion::normalize()
{
	float epsilon = 1e-7;
  float n = magnitude() + epsilon;
  c[0] /= n;
  c[1] /= n;
  c[2] /= n;
  c[3] /= n;
}

float Quaternion::magnitude() const
{ 
	float epsilon = 1e-7;
	return sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2] + c[3]*c[3] + epsilon);
}

Quaternion Quaternion::operator * (Quaternion other)
{
	Quaternion temp;
	
  temp.c[0] = c[0] * other.c[0] - c[1] * other.c[1] - c[2] * other.c[2] - c[3] * other.c[3];
  temp.c[1] = c[0] * other.c[1] + c[1] * other.c[0] + c[2] * other.c[3] - c[3] * other.c[2];
  temp.c[2] = c[0] * other.c[2] - c[1] * other.c[3] + c[2] * other.c[0] + c[3] * other.c[1];
  temp.c[3] = c[0] * other.c[3] + c[1] * other.c[2] - c[2] * other.c[1] + c[3] * other.c[0];

  return temp;
}

Quaternion Quaternion::operator * (float scaler)
{
	Quaternion temp;
	
	temp.c[0] = c[0] * scaler;
	temp.c[1] = c[1] * scaler;
	temp.c[2] = c[2] * scaler;
	temp.c[3] = c[3] * scaler;
	
	return temp;
}

Quaternion Quaternion::operator + (Quaternion other)
{
	Quaternion temp;
	
	temp.c[0] = c[0] + other.c[0];
	temp.c[1] = c[1] + other.c[1];
	temp.c[2] = c[2] + other.c[2];
	temp.c[3] = c[3] + other.c[3];
	
	return temp;
}

Quaternion Quaternion::operator - (Quaternion other)
{
	Quaternion temp;
	
	temp.c[0] = c[0] - other.c[0];
	temp.c[1] = c[1] - other.c[1];
	temp.c[2] = c[2] - other.c[2];
	temp.c[3] = c[3] - other.c[3];
	
	return temp;
}
