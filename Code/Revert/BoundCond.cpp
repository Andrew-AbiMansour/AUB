#include "Common.h"
#include "Matrix.h"

typedef Vector<double> vector; 
typedef Matrix<double> matrix;
typedef Vector<Matrix<double> > tensor;

void FixRHS(vector& RHS) {
	int nBC; double uBC;
	ifstream file;
	file.open("Dirichlet.dat");
	file >> nBC >> uBC;
	
	Vector<int> BoundPts(nBC);
	
	for(int i = 0; i < nBC; i++)
		file >> BoundPts(i);
		
	for(int i = 0; i < nBC; i++)
		RHS(BoundPts(i)) = uBC; 
}
