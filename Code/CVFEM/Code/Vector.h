#include "Main.h"
#pragma once

#ifndef VECTOR_H
#define VECTOR_H

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
	Vector<Type> operator()(const int i, const int j);
	template<typename T>
	friend ostream& operator<<(ostream& out, Vector<T>& p);
	Vector<Type>& operator()(Vector<Type>& p ,const int i, const int j);

	void Write(char* fname);

	/*** Scalar Ops ***/
	Vector<Type>& operator*=(Type scalar);
	Vector<Type>& operator+=(Type scalar);
	Vector<Type>& operator-=(Type scalar);
	Vector<Type>& operator/=(Type scalar);
	Vector<Type> operator-();
	Vector<Type> operator+();
	Vector<Type> test();
	template<typename T>
	friend Vector<T> pow(Vector<T>& p, double n);
	Type operator&&(Vector<Type>);
	Type Trapezoidal(Vector<Type>&);

	Type Amax();
	Type Amin();
	Type Amean();
	Type max();
	Type min();
	Type Norm();
	int Imax() ;

	template<typename T>
	friend Vector<T> operator*(T, Vector<T>);
	template<typename T>
	friend Vector<T> operator+(T, Vector<T>);
	template<typename T>
	friend Vector<T> operator-(T, Vector<T>);
	template<typename T>
	friend Vector<T> operator/(T, Vector<T>);

	template<typename T>
	friend Vector<T> operator+(Vector<T>,T);
	template<typename T>
	friend Vector<T> operator-(Vector<T>, T);
	template<typename T>
	friend Vector<T> operator*(Vector<T>, T);
	template<typename T>
	friend Vector<T> operator/(Vector<T>, T);

	/*** Vector Ops ***/

	template<typename T>
	friend Vector<T> operator+(Vector<T>, Vector<T>);
	template<typename T>
	friend Vector<T> operator-(Vector<T>, Vector<T>);
	template<typename T>
	friend Vector<T> operator*(Vector<T>, Vector<T>);
	template<typename T>
	friend Vector<T> operator/(Vector<T>, Vector<T>);
	
	Vector<Type> diff();
	Vector<Type>& operator+=(Vector<Type>);
	Vector<Type>& operator-=(Vector<Type>);
	Vector<Type>& operator*=(Vector<Type>);
	Vector<Type>& operator/=(Vector<Type>);

	Vector<Type>& operator=(Vector<Type>);
	Vector<Type>& operator=(Type);
	~Vector() { delete[] vec; }
};

#endif



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
	Type max = - pow(10.0,20.0);

	for(int i = 0; i < length; i++)
		max = max < fabs(vec[i]) ? fabs(vec[i]) : max;

	return max;
}

template<typename Type>
Type Vector<Type>::max() {
	Type max = - pow(10.0,20.0);

	for(int i = 0; i < length; i++)
		max = max < vec[i] ? vec[i] : max;

	return max;
}

template<typename Type>
Type Vector<Type>::min() {
	Type min = pow(10.0,20.0);

	for(int i = 0; i < length; i++)
		min = min > vec[i] ? vec[i] : min;

	return min;
}

template<typename Type>
Type Vector<Type>::Amin() {
	Type min = pow(10.0,20.0);

	for(int i = 0; i < length; i++)
		min = min > fabs(vec[i]) ? fabs(vec[i]) : min;

	return min;
}

template<typename Type>
Type Vector<Type>::Amean() {
	Type mean = 0.0;

	for(int i = 0; i < length; i++)
		mean += fabs(vec[i]);

	return mean/length;
}

template<typename Type>
int Vector<Type>::Imax() {
	Type max = 0.0;
	int j = 0;

	for(int i = 0; i < length; i++) {
		j = max < fabs(vec[i]) ? i : j;
		max = max < fabs(vec[i]) ? fabs(vec[i]) : max;
	}

	return j;
}

template<typename Type>
Type Vector<Type>::Norm() {
	Type norm = 0.0;

	for(int i = 0; i < length; i++)
		norm += pow(vec[i],2.0);

	return sqrt(norm);
}

					/// Friend Operators ///

template<typename T>
Vector<T> operator+(T scalar, Vector<T> p) {
	
	T* temp = new T[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) + scalar;

	return Vector<T>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator+(Vector<Type> p, Type scalar) {
	
	Type *temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) + scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator/(Type scalar, Vector<Type> p) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = scalar / p(i);

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator/(Vector<Type> p, Type scalar) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) / scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator-(Type scalar, Vector<Type> p) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = - p(i) + scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator-(Vector<Type> p, Type scalar) {
	Type *temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) - scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator*(Type scalar, Vector<Type> p) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) * scalar;

	return Vector<Type>(p.GetLength(),temp);
}

template<typename Type>
Vector<Type> operator*(Vector<Type> p, Type scalar) {
	
	Type* temp = new Type[p.GetLength()];
	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) * scalar;

	return Vector<Type>(p.GetLength(),temp);
}

				/// Dot Product ///

template<typename Type>
Type Vector<Type>::operator&&(Vector<Type> p) {
	assert(p.GetLength() == length);
	Type output = 0;

	for(int i = 0; i < length; i++)
		output += p(i) * vec[i];

	return output;
}
	

			/// Integration Int ///

template<typename Type>
Type Vector<Type>::Trapezoidal(Vector<Type> &p) {
	assert(p.GetLength() == length);
	Type output = 0;

	for(int i = 0; i < length - 1; i++)
		output += (p(i) + p(i+1)) * vec[i];

	return 0.5 * output;
}


/*** Vectors Ops ***/
template<typename Type>
Vector<Type> operator+(Vector<Type> p, Vector<Type> u) {
	assert(p.GetLength() == u.GetLength());

	Type* temp = new Type[p.GetLength()];

	for(int i = 0; i < p.GetLength(); i++)
		temp[i] = p(i) + u(i);

	return Vector<Type>(p.GetLength(),temp);
}
	
template<typename Type>
Vector<Type> operator-(Vector<Type> p, Vector<Type> u) {
	assert(p.GetLength() == u.GetLength());

	Type* temp = new Type[p.length];
	for(int i = 0; i < p.length; i++)
		temp[i] = p(i) - u(i);

	return Vector<Type>(p.length,temp);
}

template<typename Type>
Vector<Type> operator*(Vector<Type> p, Vector<Type> u) {
	int size = p.GetLength();
	assert(size == u.GetLength());

	Type* temp = new Type[size];
	for(int i = 0; i < size; i++)
		temp[i] = p(i) * u(i);

	return Vector<Type>(size,temp);
}

template<typename Type>
Vector<Type> operator/(Vector<Type> p, Vector<Type> u) {
	int size = p.GetLength();
	assert(size == u.GetLength());

	Type* temp = new Type[size];
	for(int i = 0; i < size; i++)
		temp[i] = p(i) / u(i);

	return Vector<Type>(size,temp);
}

template<typename Type>
Vector<Type>& Vector<Type>::operator+=(Vector<Type> p) {
	assert(p.GetLength() == length);

	for(int i = 0; i < length; i++)
		vec[i] += p(i);

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator-=(Vector<Type> p) {
	assert(p.GetLength() == length);

	for(int i = 0; i < length; i++)
		vec[i] -= p(i);

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator*=(Vector<Type> p) {
	assert(p.GetLength() == length);

	for(int i = 0; i < length; i++)
		vec[i] *= p(i);

	return *this;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator/=(Vector<Type> p) {
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
Vector<Type>::Vector(const Vector<Type> &p) : length(p.GetLength()), 
					 vec(p.GetLength() ? new Type[p.GetLength()] : 0) {
	assert(&p != this);
	for(int i = 0; i < p.GetLength(); i++)
		vec[i] = p(i);
}

template<typename Type>
Vector<Type>::Vector(int dim, Type* p) : length(dim) {
	vec = p;
}

/*** ***/

template<typename Type>
Vector<Type>& Vector<Type>::operator=(Vector<Type> p) {
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
Vector<Type> Vector<Type>::operator()(const int i, const int j) {
	const int N_temp = j - i + 1;
	Vector<Type> p(N_temp,.0);

	for(int k = 0; k < N_temp; k++)
		p(k) = vec[i+k];

	return p;
}

template<typename Type>
Vector<Type>& Vector<Type>::operator()(Vector<Type>& p ,const int i, const int j) {
	int count = 0;
	for(int k = i; k <= j; k++) {
		vec[k] = p(count);
		count++;
	}

	return *this;
}

template<typename Type>
ostream& operator<<(ostream& out, Vector<Type>& p) {
	
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
Vector<Type> pow(Vector<Type>& p, double n) {
	Vector<Type> v(p.GetLength());

	for(int i = 0; i < p.GetLength(); i++)
		v(i) = pow(p(i),n);

	return v;
}
