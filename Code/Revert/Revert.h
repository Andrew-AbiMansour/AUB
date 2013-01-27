#include <iostream>
#include <cmath>
#include "assert.h"
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib> // for GCC 4.3 and above
#include <cstring> // for GCC 4.3 and above
#pragma once
#define pi 3.141592653589793

using namespace std;

const double Da = 0.1, Db = 0.1, Dc = 0.0001, a0 = 100.0,
	         b0 = 5.0, kr = pow(10.0,-6.0), c1 = 0.1, c2 = 0.25, c3 = 0.4,
			 alfa = 10.0, beta = 10.0, gama = 10.0, delta = 10.0, 
			 velocity1 = 0.01, velocity2 = -velocity1, Nx = 3,
			 Ny = 3001, hx = 0.1, hy = 0.1;

const int Nc = 5;
