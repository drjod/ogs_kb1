#ifndef GAUSS_H
#define GAUSS_H

//#include <stdexcept>

// use slice function
/*template<typename T>
Row<T, 1> _slice(const DVec& v, const size_type& i, const size_type&)
{
	return v.slice(i);
}*/

/*
template<typename T>
Matrix_ref<T, 1>  __slice(Matrix_ref<T, 1> v, size_type i, size_type n)
{
	return v(slice(i, n));
}
*/



void classical_elimination(DMat& A, DVec& b)
{
	const size_type n = A.dim1();

	for(size_type j = 0; j != n-1; ++j)
	{
		const double pivot = A(j, j);		
		Slice_iter<double> sj = A[j];
		std::advance(sj, j);
		if(pivot == 0) throw std::runtime_error("Pivot is zero");
		for(size_type i = j+1; i != n; ++i)
		{
			const double mult = A(i, j) / pivot;
			Slice_iter<double> si = A[i];
			std::advance(si, j);

			si = scale_and_add(sj, -mult, si);
			b[i] -= mult * b[j]; 
		}
	}
}

DVec back_substitution(const DMat& A, const DVec& b)
{
	const size_type n = A.dim1();
	DVec x(n);

	for(size_type i = n-1; i >= 0; --i)
	{
		Cslice_iter<double> cxs(&x, std::slice(0, n, 1));
		Cslice_iter<double> cAi = A[i];
		const double s = b[i] - dot_product(cAi, cxs);
		++cxs, ++cAi;	

		if(const double m = A(i, i))
			x[i] = s/m;
		else
			throw std::runtime_error("Backsubstitution failure");
	}
	return x;
}

#endif
