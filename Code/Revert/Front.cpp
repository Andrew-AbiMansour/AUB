#include "Common.h"
#include "Vector.h"
#include "DataStruct.h"

double GlobRate(const DataStructure& Domain, const Vector<double>& R) {
	double GR = .0;

	for(int i = 1; i < Domain.GetYNpts() - 2; i++)
		GR += (R(Domain(i,1)) + R(Domain(i+1,1)));

	return 0.5*hy*GR;
}
