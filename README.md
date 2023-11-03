# Madgwick-MARG
An implementation of the Madgwick Quaternion Update filter, using a MARG (Magnetic, Angular Rate, and Gravity) sensor array.

For the filter to converge after startup, a higher beta filter gain is necesary (tested with `beta = 2.5` for 10 seconds, then reduced to `beta = 0.1`).

The filter gain zeta describes the amount of gyroscope bias drift compensation (tested with `zeta = 0.0015` as the gyroscope is well calibrated).

**Examples:**

```cpp
#include "MARG.h"
```
Create a MARG object:
```cpp
MARG mymarg = MARG(beta, zeta);
```
Run filter iteration (angular rate in deg/s, mag and accel are normalized later so units don't matter):
```cpp
mymarg.updateMARG(Dt, magnetic_flux[3], angular_rates[3], acceleration[3], total_angle[3]);
```
The function doesn't return anything, but changes the values in specified pointer `total_angle[3]` to the most recent attitude estimation in degrees.

Change the beta filter gain:
```cpp
mymarg.beta = 0.1f;
```
Access the calculated quaternion of orientation:
```cpp
#include "quaternions.h"
Quaternion orientation_quaternion;
orientation_quaternion = mymarg.SEq;
```
