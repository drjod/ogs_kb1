#ifndef STRU3_MATRIX_H
#define STRU3_MATRIX_H

#define CPP11

#include <iostream>
#include <valarray>
#include <algorithm>
#include <numeric>	// for std::inner_product
#include <cassert>
#include <iterator>
#include <stdexcept>

#ifdef CPP11
	#include <initializer_list>
#endif


namespace stru3
{


template<typename T> class Matrix;
typedef Matrix<double> DMat;
typedef std::valarray<double> DVec;

typedef long size_type;


#include "stru3_sliceIter.h"
//-------------------------------------------------------------


template<typename T>
class Matrix
{
	std::valarray<T>* v;	// stores elements by column as described in 22.4.5
	size_type d1, d2;	// d1 == number of rows, d2 == number of columns
public:
	Matrix(const size_type& x, const size_type& y);		// note: no default constructor
	Matrix(const Matrix&);
	Matrix& operator=(const Matrix&);
#ifdef CPP11
	Matrix(Matrix&&);
	Matrix& operator=(Matrix&&);
	//Matrix(const std::initializer_list<std::initializer_list<T> >& nested_list);
#endif
	Matrix& operator=(const T& d) { *v = d; return *this; }
	~Matrix();
	
	size_type size() const { return d1*d2; }
	size_type dim1() const { return d1; }
	size_type dim2() const { return d2; }

	Slice_iter<T> row(const size_type& i);
	Cslice_iter<T> row(const size_type& i) const;

	Slice_iter<T> column(const size_type& i);
	Cslice_iter<T> column(const size_type& i) const;

	T& operator()(const size_type& x, const size_type& y);					// Fortran-style subscripts
	T operator()(const size_type& x, const size_type& y) const;

	Slice_iter<T> operator()(const size_type& i) { return row(i); }
	Cslice_iter<T> operator()(const size_type& i) const { return row(i); }

	Slice_iter<T> operator[](const size_type& i) { return row(i); }	// C-style subscript
	Cslice_iter<T> operator[](const size_type& i) const { return row(i); }

	Matrix& operator*=(T);

	std::valarray<T>& array() { return *v; }
};

template<typename T>
inline Slice_iter<T> Matrix<T>::row(const size_type& i)
{
	return Slice_iter<T>(v,std::slice(i*d2,d2,1));
}

template<typename T>
inline Cslice_iter<T> Matrix<T>::row(const size_type& i) const
{
	return Cslice_iter<T>(v,std::slice(i*d2,d2,1));
}

template<typename T>
inline Slice_iter<T> Matrix<T>::column(const size_type& i)
{
	return Slice_iter<T>(v,std::slice(i,d1,d2));
}

template<typename T>
inline Cslice_iter<T> Matrix<T>::column(const size_type& i) const
{
	return Cslice_iter<T>(v,std::slice(i,d1,d2));
}

template<typename T>
Matrix<T>::Matrix(const size_type& x, const size_type& y)
{
	d1 = x;
	d2 = y;
	v = new std::valarray<T>(x*y);
}

template<typename T>
Matrix<T>::Matrix(const Matrix& m)
{
	v = new std::valarray<T>(m.dim1()*m.dim2());
	*v = *m.v;
	d1 = m.dim1();
	d2 = m.dim2();
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
	if(&m != this)
	{
		delete v;
		v = new std::valarray<T>(m.dim1()*m.dim2());
		*v = *m.v;
		d1 = m.dim1();
		d2 = m.dim2();
	}
	return *this;
}

#ifdef CPP11
template<typename T>
Matrix<T>::Matrix(Matrix&& m)
{
	*v = std::move(*m.v);
	d1 = m.dim1();
	d2 = m.dim2();
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>&& m)
{
	if(&m != this)
	{
		delete v;
		*v = std::move(*m.v);
		d1 = m.dim1();
		d2 = m.dim2();
	}
	return *this;
}

/*
template<typename T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T> >& nested_list)
{
	d1 = nested_list.size();
	d2 = nested_list.begin()->size();
	//const T* elems = new T[d1*d2];
	//elems = nested_list.begin()->begin();
	
	//v = new std::valarray<T>(elems, d1*d2);
	size_type i = 0;
	v = new std::valarray<T>(d1*d2);
	for(auto list: nested_list)
		for(auto elem: list)
			(*v)[i++] = elem;
}

*/
#endif

template<typename T>
Matrix<T>::~Matrix()
{
	delete v;
}

template<typename T>
T& Matrix<T>::operator()(const size_type& x, const size_type& y)
{
	return row(x)[y];
}

template<typename T>
T Matrix<T>::operator()(const size_type& x, const size_type& y) const
{
	return row(x)[y];
}


//-------------------------------------------------------------

template<typename T>
Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b)
{
	Matrix<T> res(a.dim1(), a.dim2());

	for(size_type x=0; x<a.dim1(); x++)
		for(size_type y=0; y<a.dim2(); y++)
			res(x, y) = a(x, y) + b(x, y);
	return res;
}	

template<typename T>
Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b)
{
	Matrix<T> res(a.dim1(), a.dim2());

	for(size_type x=0; x<a.dim1(); x++)
		for(size_type y=0; y<a.dim2(); y++)
			res(x, y) = a(x, y) - b(x, y);
	return res;
}	

//-------------------------------------------------------------




template<typename T>
T mul(const Cslice_iter<T>& v1, const std::valarray<T>& v2)
{
	T res = 0;
	for (size_type i = 0; i<v2.size(); i++) res+= v1[i]*v2[i];
	return res;
}


template<typename T>
std::valarray<T> operator*(const Matrix<T>& m, const std::valarray<T>& v)
{
	if (m.dim2()!=v.size()) std::cerr << "wrong number of elements in m*v\n";

	std::valarray<T> res(m.dim1());
	for (size_type i = 0; i<m.dim1(); i++) res[i] = mul(m.row(i),v);
	return res;
}


// alternative definition of m*v

//std::valarray<T> operator*(const Matrix& m, std::valarray<T>& v)
template<typename T>
std::valarray<T> mul_mv(const Matrix<T>& m, std::valarray<T>& v)
{
	if (m.dim2()!=v.size()) std::cerr << "wrong number of elements in m*v\n";

	std::valarray<T> res(m.dim1());

	for (size_type i = 0; i<m.dim1(); i++) {
		const Cslice_iter<T>& ri = m.row(i);
		res[i] = std::inner_product(ri,ri.end(),&v[0],T(0));
	}
	return res;
}



template<typename T>
std::valarray<T> operator*(std::valarray<T>& v, const Matrix<T>& m)
{
	if (v.size()!=m.dim1()) std::cerr << "wrong number of elements in v*m\n";

	std::valarray<T> res(m.dim2());

	for (size_type i = 0; i<m.dim2(); i++) {
		const Cslice_iter<T>& ci = m.column(i);
		res[i] = std::inner_product(ci,ci.end(),&v[0],T(0));
	}
	return res;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(T d)
{
	(*v) *= d;
	return *this;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m)
{
	os << '(';
	for(size_type x=0; x<m.dim1(); x++)
	{
		os << '(';
		for(size_type y=0; y<m.dim2(); y++)
		{		
			os<<m[x][y];
			if(y != m.dim2() -1) os << " ";
		}
		os << ')';
	}
	os << ')';
	return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::valarray<T>& v)
{
	os << '(';
	for(size_type x=0; x<v.size(); ++x)
	{
		os << v[x];
		if(x != v.size() -1) os << " ";
	}
	os << ')';
	return os;
}

//#include "stru3_gauss.h"

}

#endif
