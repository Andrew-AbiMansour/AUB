
#ifndef SPARSE_H
#define SPARSE_H

#include "cs.h"
#include "Vector.h"
#include "assert.h"

template<typename T> class Sparse {
	int* rows;
	int* cols;
	double* vals;
	int nnz, nnz_red, n_rows, n_cols;
	cs* Operator;

public:
	Sparse(int N, int M, int NNZ) : nnz(NNZ), n_rows(N), n_cols(M), 
		rows(new int[NNZ]), cols(new int[NNZ]), vals(new double[NNZ]) {}
	Sparse(int N, int NNZR) : nnz(NNZR), n_rows(N), n_cols(N), 
		rows(new int[NNZR]), cols(new int[N]), vals(new double[NNZR]) {}
	Sparse(int, double);
	Sparse(const Sparse<T>& S);
	Sparse(Vector<T>& v);
	Sparse(Vector<int> rows, Vector<int> cols, Vector<T> vals, int, int);

	~Sparse() { delete[] rows; delete[] cols; delete[] vals; cs_spfree(Operator); }

	int& NNZ() { return nnz; }
	int& NNZ_R() { return nnz_red; }
	int& ROWS() { return n_rows; }
	int& COLS() { return n_cols; }

	T& Val(const int i) { return vals[i]; }
	int& Row(const int i) { return rows[i]; }
	int& Col(const int i) { return cols[i]; }
	void CreateSparse();
	void CreateJacobian();
	Sparse<T> transpose();
	cs* ReturnOp() const { return Operator; }

	Sparse<T> operator+(const T);
	Sparse<T> operator-(const T);
	Sparse<T> operator*(const T);
	Sparse<T> operator/(const T);

	T& operator()(const int i, const int j);
	Sparse<T>& operator&&(Vector<T>& v);

	Sparse<T> operator+();
	Sparse<T> operator-();
	Sparse<T>& Diag(T &);
	Sparse<T>& Diag(Vector<T> &);
	Vector<T> Diag();

	template<typename Type>
	friend Sparse<Type> operator+(const Sparse<Type>&, const Sparse<Type>&);
	template<typename Type>
	friend Sparse<Type> operator-(const Sparse<Type>&, const Sparse<Type>&);
	template<typename Type>
	friend Sparse<Type> operator*(Sparse<Type>, Sparse<Type>);
	Vector<T> operator*(Vector<T>);
	Vector<T> operator/(const Vector<T>&);
	Vector<T> PCCG(const Vector<T>&, Vector<T>);

	template<typename Type>
	friend ostream& operator<<(ostream&, Sparse<T>&);
	Sparse<T>& operator=(Sparse<T>); 
	Vector<T> LoadJacobi();
	Sparse<T> LoadICChol();
	Vector<T> BiCG(const Vector<T>& RHS);
};

template<typename T>
Sparse<T>::Sparse(int N = 0, double num = 0.0) : nnz(N), n_rows(N), n_cols(N), 
		rows(new int[N]), cols(new int[N]), vals(new double[N]) {

	Operator = cs_spalloc(N,N,N,1,1);
	for(int i = 0; i < N; i++)
		cs_entry(Operator,i,i,num);

	Operator = cs_compress(Operator);

	nnz_red = Operator->nzmax;
}

template<typename T>
Sparse<T>::Sparse(Vector<int> rows, Vector<int> cols, Vector<T> vals, int N, int M) : nnz(N*M), n_rows(N), n_cols(M), 
		rows(new int[rows.GetLength()]), cols(new int[cols.GetLength()]), vals(new double[vals.GetLength()]) {

	assert(vals.GetLength() == rows.GetLength());
	assert(cols.GetLength() == vals.GetLength());

	Operator = cs_spalloc(N,M,N*M,1,1);
	for(int i = 0; i < N*M; i++)
		cs_entry(Operator,rows(i),cols(i),vals(i));

	Operator = cs_compress(Operator);

	nnz_red = Operator->nzmax;
}

template<typename T>
Sparse<T>::Sparse(const Sparse<T>& S) : nnz(S.NNZ_R()), n_rows(S.ROWS()), n_cols(S.COLS()), 
rows(new int[S.NNZ_R()]), cols(new int[S.NNZ_R()]), vals(new T[S.NNZ_R()]) {

	int count = 0;

	for(int e = 0; e < n_rows; e++) {
		int entries = S.ReturnOp()->p[e+1];

		while(count < entries) {
			rows[count] = S.ReturnOp()->i[count];
			cols[count] = e;
			vals[count] = S.ReturnOp()->x[count];
			count++;
		}
	}

	CreateSparse();
}

template<typename T>
Sparse<T>::Sparse(Vector<T>& v) : nnz(v.GetLength()), n_rows(v.GetLength()), n_cols(v.GetLength()), 
rows(new int[v.GetLength()]), cols(new int[v.GetLength()]), vals(new T[v.GetLength()]) {

	int count = 0;

	for(int e = 0; e < n_rows; e++) {

			rows[e] = e;
			cols[e] = e;
			vals[e] = v(e);
	}

	CreateSparse();
}


template<typename T>
void Sparse<T>::CreateSparse() {

	Operator = cs_spalloc(n_cols,n_rows,nnz,0,1);
	free(Operator->i); free(Operator->p);  free(Operator->x);

	Operator->i = rows;
	Operator->p = cols;
	Operator->x = vals;
	Operator->nz = nnz;

	Operator = cs_compress(Operator);
	cs_dupl(Operator);

	nnz_red = Operator->nzmax;
}

template<typename T>
void Sparse<T>::CreateJacobian() {

	Operator = cs_spalloc(n_cols,n_rows,nnz,1,1);
	free(Operator->i); free(Operator->p);  free(Operator->x);

	Operator->i = rows;
	Operator->p = cols;
	Operator->x = vals;
	Operator->nz = nnz;

	nnz_red = Operator->nzmax;
}

template<typename T>
Sparse<T>& Sparse<T>::operator=(Sparse<T> S) {

	//Sparse::~Sparse();
	delete[] rows; delete[] cols; delete[] vals; cs_spfree(this->ReturnOp());

	nnz = S.NNZ(); 
	n_rows = S.ROWS(); 
	n_cols = S.COLS();

	rows = new int[nnz];
	cols = new int[nnz];
	vals = new T[nnz];

	int count = 0;

	for(int e = 0; e < n_rows; e++) {
		int entries = S.ReturnOp()->p[e+1];

		while(count < entries) {
			rows[count] = S.ReturnOp()->i[count];
			cols[count] = e;
			vals[count] = S.ReturnOp()->x[count];
			count++;
		}
	}

	CreateSparse();
	return *this;
}

template<typename T>
Sparse<T> Sparse<T>::operator+(const T a) {
	Sparse<T> C = *this;

	for(int i = 0; i < Operator->nzmax; i++)
		 C.ReturnOp()->x[i] = Operator->x[i] + a;

	return C;
}

template<typename T>
Sparse<T> Sparse<T>::operator-(const T a) {
	Sparse<T> C = *this;

	for(int i = 0; i < Operator->nzmax; i++)
		 C.ReturnOp()->x[i] = Operator->x[i] - a;

	return C;
}

template<typename T>
Sparse<T> Sparse<T>::operator*(const T a) {
	Sparse<T> C = *this;

	for(int i = 0; i < Operator->nzmax; i++)
		 C.ReturnOp()->x[i] = Operator->x[i] * a;

	return C;
}

template<typename T>
Sparse<T> Sparse<T>::operator/(const T a) {
	Sparse<T> C = *this;

	for(int i = 0; i < Operator->nzmax; i++)
		 C.ReturnOp()->x[i] = Operator->x[i] / a;

	return C;
}

template<typename T>
Sparse<T> operator+(const Sparse<T>& M, const Sparse<T>& N) {
	Sparse<T> C(M.n_rows,M.n_cols,M.nnz);
	C.Operator = cs_add(M.ReturnOp(),N.ReturnOp(),1.0,1.0);

	C.NNZ_R() = C.ReturnOp()->nzmax;

	return C;
}

template<typename T>
Sparse<T> operator-(const Sparse<T>& M, const Sparse<T>& N) {
	Sparse<T> C(M.n_rows,M.n_cols,M.nnz);
	C.Operator = cs_add(M.ReturnOp(),N.ReturnOp(),1.0,-1.0);

	C.NNZ_R() = C.ReturnOp()->nzmax;

	return C;
}

template<typename T>
Sparse<T> operator*(Sparse<T> M, Sparse<T> N) {
	assert(M.n_cols == N.n_rows);
	cs* Operator = cs_multiply(M.ReturnOp(),N.ReturnOp());

	Sparse<T> C(M.ROWS(),N.COLS(),Operator->nzmax);
	C.Operator = Operator;

	C.NNZ_R() = Operator->nzmax;
	return C;
}

template<typename T>
ostream& operator<<(ostream& out, Sparse<T>& S) {
	cs_print(S.ReturnOp(),0);
	return out;
}

template<typename T>
Vector<T> Sparse<T>::operator*(Vector<T> v) {
	Vector<T> p(v.GetLength());
	cs_gaxpy(Operator,v.GetAddress(),p.GetAddress());

	return p;
}

template<typename T>
Vector<T> Sparse<T>::operator/(const Vector<T>& v) {
	cs_lusol(1,Operator,v.GetAddress(),0);
	//cs_cholsol(0,Operator,v.GetAddress());
	Vector<double> p(v);

	return p;
}

template<typename T>
Sparse<T> Sparse<T>::transpose() {
	Sparse<T> trans(n_cols,n_rows,Operator->nzmax);
	trans.Operator = cs_transpose(Operator,1);
	trans.NNZ_R() = Operator->nzmax;

	return trans;
}

template<typename T>
Vector<T> Sparse<T>::PCCG(const Vector<T>& b, Vector<T> x0) {

	Vector<T> r0 = b - (*this) * x0;
	if(r0.Norm() == .0) {
		Vector<T> z(x0.GetLength(),pow(10.0,-7.0));
		r0 = b - (*this) * z;
		cout << "Initial residual in CG iz zero!" << endl;
	}

	//sparse PC = LoadICChol();
	Vector<T> PC = LoadJacobi();
	Vector<T> z0 = PC * r0;
	Vector<T> p0 = z0;
	T error = 1.0, tol = pow(10.0,-6.0);

	int iters = 0;
	while( error > tol ) {
		Vector<T> Ap0 = (*this) * p0;
		T alfa = (r0 && z0) / (p0 && Ap0 );
		x0 += alfa * p0;
		T rtr = (z0 && r0);
		r0 -= alfa * Ap0;
		z0 = PC * r0;

		error = sqrt(r0 && r0);
		T beta = (z0 && r0) / rtr;
		p0 = beta * p0 + z0; 
		iters++;
	}
	//cout << "CG converged in " << iters << endl;
	return x0;
}  

template<typename T>
Sparse<T> Sparse<T>::LoadICChol() {

	css *order = cs_schol(1,Operator);
	csn* Chol = cs_chol(Operator,order);
	cs* IChol = cs_multiply(Chol->L,Chol->U); 

	cout << IChol->nzmax << endl;
	Sparse<T> C(n_rows,n_cols,IChol->nzmax);
	C.Operator = IChol;

	C.NNZ_R() = IChol->nzmax;

	return C;
}

template<typename T>
Vector<T> Sparse<T>::LoadJacobi() {
	int index1 = 0, index2 = 0;
	Vector<T> PC(n_cols);
	
	for(int k = 0; k < n_cols ; k++) {

		index2 = Operator->p[k+1];

		for(int j = index1; j < index2; j++)

			if( Operator->i[j] == k )
				PC(k) = 1.0 / Operator->x[j];

		index1 = index2;
	}
	return PC;
}

template<typename T>
Sparse<T> Sparse<T>::operator+() {
	return *this;
}

template<typename T>
Sparse<T> Sparse<T>::operator-() {
	Sparse<T> S = *this;

	for(int i = 0; i < Operator->nzmax; i++)
		 S.ReturnOp()->x[i] = - Operator->x[i];

	return S;
}

template<typename T>
Vector<T> Sparse<T>::BiCG(const Vector<T>& RHS)
{
	int NVERTS = RHS.GetLength();
	Vector<T> r0(NVERTS), r(NVERTS), p(NVERTS), v(NVERTS), t(NVERTS), s(NVERTS), x(NVERTS);
	Vector<T> PC = LoadJacobi();
	r0 = RHS;
	r = r0;
	int iters = 0;

	double rho0 = 1.0, rho1, w0 = 1.0, alfa = 1.0, beta, Norm_r = 1.0, tol = pow(10.0,-10);

	while( Norm_r > tol ) {
		iters++;
		rho1 = r0 && r;

		beta = (rho1 / rho0) * (alfa / w0);

		p = r + beta * (p - w0 * v);

		Vector<T> y = PC * p;

		v = (*this) * y;

		alfa = rho1 / (r0 && v);
		s = r - alfa * v;

		Vector<T> z = PC * s;

		t = (*this) * z;
		Vector<T> temp = PC * t;

		w0 = (temp && (PC * s)) / (temp && temp);

		x += alfa * y + w0 * z;
		r = s - w0 * t;

		rho0 = rho1;

		Norm_r = r.Norm();
	}
	return x;
}

template<typename T>
Sparse<T>& Sparse<T>::Diag(T& s) {
	int index1 = 0, index2 = 0;

	for(int k = 0; k < n_cols; k++) {
		index2 = Operator->p[k+1];
		for(int j = index1; j < index2; j++)
			if(Operator->i[j] == k)
				Operator->x[j] += s;
		index1 = index2;
	}

	return *this;
}

template<typename T>
Sparse<T>& Sparse<T>::Diag(Vector<T>& s) {
	int index1 = 0, index2 = 0;

	for(int k = 0; k < n_cols; k++) {
		index2 = Operator->p[k+1];
		for(int j = index1; j < index2; j++)
			if(Operator->i[j] == k)
				Operator->x[j] += s(k);
		index1 = index2;
	}

	return *this;
}

template<typename T>
Vector<T> Sparse<T>::Diag() {
	int index1 = 0, index2 = 0;
	assert(n_rows == n_cols);

	Vector<T> V(n_rows,.0);

	for(int k = 0; k < n_cols; k++) {
		index2 = Operator->p[k+1];
		for(int j = index1; j < index2; j++)
			if(Operator->i[j] == k)
				V(k) += Operator->x[j];
		index1 = index2;
	}

	return V;
}

template<typename T>
T& Sparse<T>::operator()(const int row, const int col) {
	int entries = Operator->p[col+1] - Operator->p[col];
	T a;

	for(int j = 0; j < entries; j++) {
		int index = Operator->p[col]+j;
	if(Operator->i[index] == row)
		return Operator->x[row]; }

	return a;
}

template<typename T>
Sparse<T>& Sparse<T>::operator&&(Vector<T>& v) {
	
	int index1 = 0, index2;
	for(int j = 0; j < v.GetLength(); j++) {
		index2 = Operator->p[j+1];
		for(int k = index1; k < index2; k++)
			Operator->x[k] *= v(j);
		index1 = index2;
	}

	return *this;
}
#endif
