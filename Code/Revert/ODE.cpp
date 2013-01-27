#include "Common.h"
#include "Matrix.h"
#include "Sparse.h"

typedef Vector<double> vector;
typedef Matrix<double> matrix;
typedef Sparse<double> sparse;

/*vector ConstructEqs(const vector& x, const matrix& RHS, const sparse& L);

sparse ConstructJacobian(const vector& x, const matrix& RHS, const sparse& Diff)
{ 
	unsigned N = x.GetLength(), nnz = 0, NNZMAX = 3*N*Nc;
	int *rows = new int[NNZMAX], 
	    *cols = new int[NNZMAX];
	    
	vector Derivative(N);
	double h = pow(10.0,-6.0), *vals = new double[NNZMAX];
	
	for(int i = 0; i < N; i++) {
		vector Eqs = ConstructEqs(x,RHS,Diff);
		x(i) += h;
		vector NewEqs = ConstructEqs(x,RHS,Diff);
		x(i) -= h;
	
		Derivative = (NewEqs - Eqs);
		
		for(int k = 0; k < N; k++)
			if(Derivative(k) != .0) {
				rows[nnz] = k;
				cols[nnz] = i;
				vals[nnz] = Derivative(k) / h;
				nnz++;
		}
	}

	assert(nnz <= NNZMAX);
	sparse Jacobian(N,N,nnz);
	
	for(int k = 0; k < nnz; k++) {
			Jacobian.Row(k) = rows[k];
			Jacobian.Col(k) = cols[k];
			Jacobian.Val(k) = vals[k];
	}
		
	Jacobian.CreateSparse();
	delete[] rows; delete[] cols; delete[] vals;
	return Jacobian;
}
*/

vector ComputeRecursiveDerivative(int order, int Gorder, int n, const vector& dt, const vector& x, vector& der) {

	double coef = dt(0), sum = dt(0);
	
	for(int i = 1; i < Gorder; i++)
		sum += dt(i);
	
	coef /= sum;
	
	if(order == 1)
		der = coef * (x(n-1) - x(n)) / (dt(n-1)); 
	else
		der =  (ComputeRecursiveDerivative(order-1,Gorder,n-1,dt,x,der) - ComputeRecursiveDerivative(order-1,Gorder,n,dt,x,der));
		
	//cout << "der : " << order << " - " << der << endl;  
	return der;
}

matrix ConstructSequence(const vector& dt, const vector& x) {
	int Order = dt.GetLength();
	matrix Sequence(Order,1);
	
	for(int i = Order; i > 0; i--)
			Sequence[i-1] = ComputeRecursiveDerivative(i,i,i,dt,x,Sequence[i-1]);
	
	return Sequence;
} 	

vector ConstructBDF(int Order, const vector& dt) {
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


vector ConstructBDFConst(int length) {
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
