#ifndef VECTOR_H
#define VECTOR_H
#include "Common.h"

template <typename Type>
class Vector
{
protected:
	Type *vec;
	int length;

public:
	Type* GetAddress() { return vec; }
	void Allocate(int, Type);
	/*** Constructors ***/
	Vector(int, Type a = 0.0);
	Vector(const Vector &);
	Vector(int, Type*);
	Vector() {}

	/*** Access Ops ***/
	int GetLength() const { return length; }
	Type& operator()(int i) const;
	template<typename T>
	friend ostream& operator<<(ostream& out, const Vector<T>& p);
	Vector<Type>& operator()(const Vector<Type>& p ,const int i, const int j);
	Vector<Type> operator()(const int i, const int j,const Vector<Type>& p);
	Vector<Type> operator()(const int i, const int j);
	
	void Write(char* fname);
	void Read(char* fname);
	
	/*** Scalar Ops ***/
	Vector<Type>& operator*=(Type scalar);
	Vector<Type>& operator+=(Type scalar);
	Vector<Type>& operator-=(Type scalar);
	Vector<Type>& operator/=(Type scalar);
	Vector<Type> operator-();
	Vector<Type> operator+();
	Vector<Type> test();
	template<typename T>
	friend Vector<T> pow(const Vector<T>& p, double n);
	Type operator&&(const Vector<Type>&);
	Type Trapezoidal(const Vector<Type>&);

	Type Amax();
	Type Amin();
	Type Amean();
	Type max();
	Type min();
	Type Norm();
	int Imax() ;
	Type Sum();

	template<typename T>
	friend Vector<T> operator*(T, const Vector<T>&);
	template<typename T>
	friend Vector<T> operator+(T, const Vector<T>&);
	template<typename T>
	friend Vector<T> operator-(T, const Vector<T>&);
	template<typename T>
	friend Vector<T> operator/(T, const Vector<T>&);

	template<typename T>
	friend Vector<T> operator+(const Vector<T>&,T);
	template<typename T>
	friend Vector<T> operator-(const Vector<T>&, T);
	template<typename T>
	friend Vector<T> operator*(const Vector<T>&, T);
	template<typename T>
	friend Vector<T> operator/(const Vector<T>&, T);

	/*** Vector Ops ***/

	template<typename T>
	friend Vector<T> operator+(const Vector<T>&, const Vector<T>&);
	template<typename T>
	friend Vector<T> operator-(const Vector<T>&, const Vector<T>&);
	template<typename T>
	friend Vector<T> operator*(const Vector<T>&, const Vector<T>&);
	template<typename T>
	friend Vector<T> operator/(const Vector<T>&, const Vector<T>&);
	
	Vector<Type> diff();
	Vector<Type>& operator+=(const Vector<Type>&);
	Vector<Type>& operator-=(const Vector<Type>&);
	Vector<Type>& operator*=(const Vector<Type>&);
	Vector<Type>& operator/=(const Vector<Type>&);

	Vector<Type>& operator=(const Vector<Type>&);
	Vector<Type>& operator=(Type);
	~Vector() { delete[] vec; }
};

/*template<typename T>
Vector<T> Vector<T>::NewtonKrylov(const Vector<T>& b, Vector<T> x0, Vector<T> (*Jdv)(const Vector<T>&, const Vector<T>&)) {

	Vector<T> r0 = b - Jdv(x0,x0);

	if(r0.Norm() == .0) {
		Vector<T> z(x0.GetLength(),pow(10.0,-7.0));
		r0 = b - (*this) * z;
		cout << "Initial residual in CG iz zero!" << endl;
	}

	Vector<T> PC = LoadJacobi();
	Vector<T> z0 = PC * r0;
	Vector<T> p0 = z0;
	T error = 1.0, tol = pow(10.0,-6.0);

	int iters = 0;
	while( error > tol ) {
		Vector<T> Ap0 = Jdv(x0,p0);
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

	return x0;
} */ 

template<typename Type>
void Vector<Type>::Allocate(int n, Type a) { 
	vec = new Type[n];
	length = n;
	(*this) = a;
}

template<typename Type>
Vector<Type> Vector<Type>::operator-() {
	Vector<Type> v = *this;

	for(int i = 0; i < length; i++)
		v(i) = - vec[i];

	return v;
}

template<typename Type>
Vector<Type> Vector<Type>::operator+() {
	Vector<Type> v = *this;
	return v;
}

/*** Scalar Ops ***/
template<typename Type>
Vector<Type>& Vector<Type>::operator+=(Type scalar) {
	
	for(int i = 0; i < length; i++)
		vec[i] += scalar;

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator/=(Type scalar) {
	
	for(int i = 0; i < length; i++)
		vec[i] /= scalar;

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator-=(Type scalar) {
	
	for(int i = 0; i < length; i++)
		vec[i] -= scalar;

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator*=(Type scalar) {
	for(int i = 0; i < length; i++)
		vec[i] *= scalar;
	return *this;
}

template<typename Type>
Type Vector<Type>::Amax() {
	Type max = vec[0];

	for(int i = 1; i < length; i++)
		max = max < fabs(vec[i]) ? fabs(vec[i]) : max;

	return max;
}

template<typename Type>
Type Vector<Type>::Sum() {
	Type sum = vec[0];

	for(int i = 1; i < length; i++)
		sum += vec[i];

	return sum;
}

template<typename Type>
Type Vector<Type>::max() {
	Type max = vec[0];

	for(int i = 1; i < length; i++)
		max = max < vec[i] ? vec[i] : max;

	return max;
}

template<typename Type>
Type Vector<Type>::min() {
	Type min = vec[0];

	for(int i = 1; i < length; i++)
		min = min > vec[i] ? vec[i] : min;

	return min;
}

template<typename Type>
Type Vector<Type>::Amin() {
	Type min = vec[0];

	for(int i = 1; i < length; i++)
		min = min > fabs(vec[i]) ? fabs(vec[i]) : min;

	return min;
}

template<typename Type>
Type Vector<Type>::Amean() {
	Type mean = fabs(vec[0]);

	for(int i = 1; i < length; i++)
		mean += fabs(vec[i]);

	return mean/length;
}

template<typename Type>
int Vector<Type>::Imax() {
	Type max = vec[0];
	int j = 0;

	for(int i = 1; i < length; i++) {
		j = max < fabs(vec[i]) ? i : j;
		max = max < fabs(vec[i]) ? fabs(vec[i]) : max;
	}

	return j;
}

template<typename Type>
Type Vector<Type>::Norm() {
	Type norm = pow(vec[0],2.0);

	for(int i = 1; i < length; i++)
		norm += pow(vec[i],2.0);

	return sqrt(norm);
}

					/// Friend Operators ///

template<typename T>
Vector<T> operator+(T scalar, const Vector<T>& p) {
	
	T* temp = new T[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) + scalar;

	return Vector<T>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator+(const Vector<Type>& p, Type scalar) {
	
	Type *temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) + scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator/(Type scalar, const Vector<Type>& p) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = scalar / p(i);

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator/(const Vector<Type>& p, Type scalar) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) / scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator-(Type scalar, const Vector<Type>& p) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = - p(i) + scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator-(const Vector<Type>& p, Type scalar) {
	Type *temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) - scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator*(Type scalar, const Vector<Type>& p) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) * scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator*(const Vector<Type>& p, Type scalar) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) * scalar;

	return Vector<Type>(p.GetLength(),temp);
}

				/// Dot Product ///

template<typename Type>
Type Vector<Type>::operator&&(const Vector<Type>& p) {
	assert(p.GetLength() == length);
	Type output = 0;

	for(int i = 0; i < length; i++)
		output += p(i) * vec[i];

	return output;
}
	

			/// Integration Int ///

template<typename Type>
Type Vector<Type>::Trapezoidal(const Vector<Type> &p) {
	assert(p.GetLength() == length);
	Type output = 0;

	for(int i = 0; i < length - 1; i++)
		output += (p(i) + p(i+1)) * vec[i];

	return 0.5 * output;
}


/*** Vectors Ops ***/
template<typename Type>
Vector<Type> operator+(const Vector<Type>& p, const Vector<Type>& u) {
	assert(p.GetLength() == u.GetLength());

	Type* temp = new Type[p.GetLength()];

	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) + u(i);

	return Vector<Type>(p.GetLength(),temp);
}
	
template<typename Type>
Vector<Type> operator-(const Vector<Type>& p, const Vector<Type>& u) {
	assert(p.GetLength() == u.GetLength());

	Type* temp = new Type[p.length];
	for(int i = 0; i < p.length; i++)
		temp[i] = p(i) - u(i);

	return Vector<Type>(p.length,temp);
}

template<typename Type>
Vector<Type> operator*(const Vector<Type>& p, const Vector<Type>& u) {
	int size = p.GetLength();
	assert(size == u.GetLength());

	Type* temp = new Type[size];
	for(int i = 0; i < size; i++)
		temp[i] = p(i) * u(i);

	return Vector<Type>(size,temp);
}

template<typename Type>
Vector<Type> operator/(const Vector<Type>& p, const Vector<Type>& u) {
	int size = p.GetLength();
	assert(size == u.GetLength());

	Type* temp = new Type[size];
	for(int i = 0; i < size; i++)
		temp[i] = p(i) / u(i);

	return Vector<Type>(size,temp);
}

template<typename Type>
Vector<Type>& Vector<Type>::operator+=(const Vector<Type>& p) {
	assert(p.GetLength() == length);

	for(int i = 0; i < length; i++)
		vec[i] += p(i);

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator-=(const Vector<Type>& p) {
	assert(p.GetLength() == length);

	for(int i = 0; i < length; i++)
		vec[i] -= p(i);

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator*=(const Vector<Type>& p) {
	assert(p.GetLength() == length);

	for(int i = 0; i < length; i++)
		vec[i] *= p(i);

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator/=(const Vector<Type>& p) {
	assert(p.GetLength() == length);

	for(int i = 0; i < length; i++)
		vec[i] /= p(i);

	return *this;
}

template<typename Type>
Vector<Type> Vector<Type>::diff() {
	Vector<Type> p(length - 1,.0);

	for(int i = 0; i < length - 1; i++)
		p(i) = vec[i+1] - vec[i];

	return p;
}

/*** Constructors ***/
template<typename Type>
Vector<Type>::Vector(int dim, Type a) : length(dim), vec(dim ? new Type[dim] : 0) { 
	for(int i = 0; i < length; i++)
		vec[i] = a;
}

template<typename Type>
Vector<Type>::Vector(const Vector<Type>& p) : length(p.GetLength()), 
					 vec(p.GetLength() ? new Type[p.GetLength()] : 0) {
	for(int i = 0; i < p.GetLength(); i++)
		vec[i] = p(i);
}

template<typename Type>
Vector<Type>::Vector(int dim, Type* p) : length(dim) {
	vec = p;
}

/*** ***/

template<typename Type>
Vector<Type>& Vector<Type>::operator=(const Vector<Type>& p) {
	assert(p.GetLength() == length);
	for(int i = 0; i < length; i++)
		vec[i] = p(i);

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator=(Type a) {
	for(int i = 0; i < length; i++)
		vec[i] = a;

	return *this;
}

/*** Access operator ***/
template<typename Type>
Type& Vector<Type>::operator()(int i) const{
	assert(i >=0 && i <length);
	return vec[i];
}

template<typename Type>
Vector<Type>& Vector<Type>::operator()(const Vector<Type>& p ,const int i, const int j) {
	assert(i >=0 && j < length);
	
	for(int k = i; k <= j; k++)
		vec[k] = p(k-i);

	return *this;
}

template<typename Type>
Vector<Type> Vector<Type>::operator()(const int i, const int j) {
	assert(i >= 0 && j < length);
	 
	Vector<Type> temp(j-i+1);
	
	for(int k = i; k <= j; k++)
		temp(k-i) = (*this)(k);
		
	return temp;
}

template<typename Type>
Vector<Type> Vector<Type>::operator()(const int i, const int j, const Vector<Type>& p) {
	assert(i >= 0 && j < p.GetLength());
	 
	Vector<Type> temp(j-i+1);
	
	for(int k = i; k <= j; k++)
		temp(k-i) = p(k);
		
	return temp;
}

template<typename Type>
ostream& operator<<(ostream& out, const Vector<Type>& p) {
	
	for(int i = 0; i < p.GetLength(); i++)
		out << p(i) << " ";

	out << endl;
	return out;
}

template<typename Type>
void Vector<Type>::Write(char* fname) {
	ofstream myfile;
  	myfile.open (fname,ios::app);

	for(int i = 0; i < length; i++)
		myfile << vec[i] << endl;

	myfile.close();
}

template<typename Type>
void Vector<Type>::Read(char* fname) {
	ifstream myfile;
  	myfile.open (fname,ios::app);

	for(int i = 0; i < length; i++)
		myfile >> vec[i] >> endl;

	myfile.close();
}

template<typename Type>
Vector<Type> pow(const Vector<Type>& p, double n) {
	Vector<Type> v(p.GetLength());

	for(int i = 0; i < p.GetLength(); i++)
		v(i) = pow(p(i),n);

	return v;
}

#endif
