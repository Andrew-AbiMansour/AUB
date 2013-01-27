#include "Vector.h"

#ifndef MATRIX_H
#define MATRIX_H

#pragma once

template<typename T> 
class Matrix: public Vector<Vector<T> > 
{
protected:
	int rows, cols;

public:
	/*** Constructors ***/
	Matrix(int dim1, int dim2, const T);
	Matrix(const Matrix<T>&);
	Matrix() {;}

	/*** Access operators ***/
	Vector<T> operator()(const int);
	Vector<T>& operator[](const int) const;
	T& operator()(const int i, const int j) const {
		return Vector<Vector<T> >::vec[i](j); }
	template<typename Type>
	friend ostream& operator<<(ostream&, Matrix<Type>&); 
	int GetRows() const { return rows; }
	int GetCols() const { return cols; }
	int GetLength() const { return Vector<T>::length; }

	Matrix<T>& operator=(T);
	Matrix<T>& operator=(Matrix<T>&);

	/*** Matrix-Matrix Operators ***/
	Matrix<T> operator+(Matrix<T>&);
	Matrix<T> operator-(Matrix<T>&);
	Matrix<T> operator*(Matrix<T>&);
	Matrix<T> operator/(Matrix<T>&);

	/*** Matrix-Vector Operators ***/
	Vector<T> operator*(Vector<T>&);
	template<typename Type>
	friend Vector<Type> operator*(Vector<Type>& p, Matrix<Type> &m);

	/*** Matrix self-operators ***/
	template<typename Type>
	Matrix<Type> trans(Matrix<Type>&);
	~Matrix() {;}
};

#endif


/*** Constructors ***/
template<typename T>
Matrix<T>::Matrix(int dim1, int dim2, const T a = 0.0) : rows(dim1), cols(dim2) {
	Vector<Vector<T> >::length = dim1;
	Vector<Vector<T> >::vec = rows ? new Vector<T>[rows] : 0;

	for(int i = 0; i < dim1; i++)
		Vector<Vector<T> >::vec[i].Allocate(dim2,a);
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>& m) {
	Vector<Vector<T> >::vec = m.GetRows() ? new Vector<T>[m.GetRows()] : 0;
	Vector<Vector<T> >::length = m.GetRows();
	//m.length = Vector<Vector<T> >::length;
	rows = Vector<Vector<T> >::length;
	cols = m.GetCols();

	for(int i = 0; i < rows; i++)
		Vector<Vector<T> >::vec[i].Allocate(cols,0.0);

	for(int i  = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			Vector<Vector<T> >::vec[i](j) = m(i,j);
}


/*** Access operators ***/
template<typename T>
Vector<T> Matrix<T>::operator()(const int i) {
	assert(i >= 0 && i < cols);

	Vector<T> temp(rows,.0);
	for(int j = 0; j < rows; j++) temp(j) = Vector<Vector<T> >::vec[j](i);

	return temp;
}

template<typename T>
Vector<T>& Matrix<T>::operator[](const int i) const {
	assert(i >= 0 && i <= rows);
	return Vector<Vector<T> >::vec[i];
}

template<typename T>
ostream& operator<<(ostream& out, Matrix<T>& m) {

	for(int i = 0; i < m.GetRows(); i++) {
		for(int j = 0; j < m.GetCols(); j++)
				out << m(i,j) << " ";
		out << endl;
	}

	return out;
}


/*** Self-Matrix operators ***/
template<typename T>
Matrix<T> trans(Matrix<T>& m) {

	Matrix<T> mt(m.GetCols(),m.GetRows(),.0);

	for(int i = 0; i < m.GetRows(); i++)
		for(int j = 0; j < m.GetCols(); j++)
			mt(j,i) = m(i,j);

	return mt;
}

/*** Matrix-Matrix operators ***/
template<typename T>
Matrix<T> Matrix<T>::operator+(Matrix<T>& m) {

	Matrix<T> temp(rows,cols);
	for(int i = 0; i < cols; i++)
		temp(i) = m(i) + Vector<Vector<T> >::vec[i];

	return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(Matrix<T>& m) {

	Matrix<T> temp(rows,cols);
	for(int i = 0; i < cols; i++)
		temp(i) = Vector<Vector<T> >::vec[i] - m(i);

	return temp;
}

/*** Matrix-vector operation(s) ***/

template<typename T>
Vector<T> Matrix<T>::operator*(Vector<T>& p) {

	Vector<T> temp(rows);
	for(int i = 0; i < cols; i++)
		temp += Vector<Vector<T> >::vec[i] * p(i);

	return temp;
}

template<typename T>
Vector<T> operator*(Vector<T>& p, Matrix<T> &m) {

	Vector<T> temp(m.GetCols());
	for(int i = 0; i < m.GetRows(); i++)
		temp += m(i) * p(i);

	return temp;
}
			
template<typename T>
Matrix<T>& Matrix<T>::operator=(T s) {
	for(int i  = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			Vector<Vector<T> >::vec[i](j) = s;

	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>& m) {
	Vector<Vector<T> >::vec = m.GetRows() ? new Vector<T>[m.GetRows()] : 0;
	Vector<Vector<T> >::length = m.GetRows();
	m.length = Vector<Vector<T> >::length;
	rows = Vector<Vector<T> >::length;
	cols = m.GetCols();

	for(int i = 0; i < rows; i++)
		Vector<Vector<T> >::vec[i].Allocate(cols,0.0);

	for(int i  = 0; i < rows; i++)
		for(int j = 0; j < cols; j++)
			Vector<Vector<T> >::vec[i](j) = m(i,j);

	return *this;
}
