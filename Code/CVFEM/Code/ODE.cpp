#include "Main.h"
#include "Vector.h"

typedef Vector<double> vector;

vector Construct_BDF(int Order, const vector& dt) {
	vector coef(Order+1,.0);
	Vector<int> occ1(Order,1), occ2(Order,0);
	coef(0) = 1.0; coef(1) = -1.0;

	for(int i = 0; i < 2; i++) {
		int num = (i%2) ? -1 : 1;

			for(int j = 1; j < Order; j++) {
				double factor = .0;

				for(int k = 0; k <= j; k++)
					factor += dt(k);
				factor = dt(0) * dt(0) / factor;

				// The 'i' is optional ... could use an if statement too //
				occ2(j) = occ1(j-1) + i * occ2(j-1);

				if( i > 0)
					coef(i) += num * (occ1(j-1) / dt(i-1) + occ2(j-1) / dt(i)) * factor;
				else
					coef(i) += num * (occ1(j-1) / dt(i)) * factor;
			}

		occ2(0) = 1;
	}

	occ1 = occ2;

	for(int i = 2; i <= Order; i++) {
		int num = (i%2) ? -1 : 1;

		occ2(i-2) = 0;

		for(int j = i - 1; j < Order; j++) {
			double factor = .0;

			for(int k = 0; k <= j; k++)
				factor += dt(k);
			factor = dt(0) * dt(0) / factor;

			occ2(j) = occ1(j-1) + occ2(j-1);
			if(i == Order)
				coef(i) += num * (occ1(j-1) / dt(i-1)) * factor;
			else
				coef(i) += num * (occ1(j-1) / dt(i-1) + occ2(j-1) / dt(i)) * factor;
		}

		occ1 = occ2;
	}

	return coef/dt(0);
}


vector Construct_BDF_const(int length) {
	vector coef(length+1,.0);
	Vector<int> occ1(length,1), occ2(length,0);
	coef(0) = 1.0; coef(1) = -1.0;

	for(int i = 0; i < 2; i++) {
		int num = (i%2) ? -1 : 1;

			for(int j = 1; j < length; j++) {
				occ2(j) = occ1(j-1) + i * occ2(j-1);
				coef(i) += num * occ2(j) / (j+1.0);
			}

		occ2(0) = 1;
	}

	occ1 = occ2;

	for(int i = 2; i <= length; i++) {
		int num = (i%2) ? -1 : 1;

		occ2(i-2) = 0;

		for(int j = i - 1; j < length; j++) {
			occ2(j) = occ1(j-1) + occ2(j-1);
			coef(i) += num * occ2(j) / (j+1.0);
		}

		occ1 = occ2;
	}

	return coef;
}
