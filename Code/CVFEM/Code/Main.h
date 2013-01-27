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

/*const double Da = pow(10.0,-4.0), Db = Da, Dc = 1.1 *  pow(10.0,-4.0), 
a0 = 1000.0,
	     b0 = 3.0, kr = pow(10.0,-5.0), c1 = 1.1, c2 = 2.0, c3 = 2.0, 
	     alfa = 0.01, beta = 0.01, gama = 0.01, delta = 0.003;

Params above lead to thick rings.
*/

/* Multi-circles params */
/*
const double Da = pow(10.0,-4.0), Db = Da, Dc = 5.0 * pow(10.0,-5.0),
     	     a0 = 500.0, b0 = 3.0, kr = pow(10.0,-5.0),
	     c1 = 1.1, c2 = 2.0, c3 = 2.0,
	     alfa = 0.01, beta = 0.01, gama = 0.01, delta = 0.003;
*/

/* Fast systems */
const double Da = 1.0, Db = Da, Dc = 1.0, a0 = 500.0,b0 = 3.0, kr = 
0.001, c1 = 1.0, c2 = 1.25, c3 = 1.5, alfa = 2.0, beta = 1.0, gama = 1.0,
 delta = 0.5;

/* Triangle
const double Da = pow(10.0,-3.0), Db = Da, Dc = 1.1 * Da, a0 = 1500,
	     b0 = 3.0, kr = 5.0 * pow(10.0,-5.0), c1 = 1.7, c2 = 2.0, 
	     c3 = 2.0, alfa = 0.1, beta = 0.1, gama = 0.1, delta = 0.1;
*/
const int Nc = 5;
