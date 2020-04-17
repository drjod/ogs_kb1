#ifndef SLICE_ITER_H
#define SLICE_ITER_H

template<typename T> class Slice_iter;
template<typename T> bool operator==(const Slice_iter<T>&, const Slice_iter<T>&);
template<typename T> bool operator!=(const Slice_iter<T>&, const Slice_iter<T>&);
template<typename T> bool operator< (const Slice_iter<T>&, const Slice_iter<T>&);


template<typename T> class Slice_iter : public std::iterator<std::forward_iterator_tag, T>
{

	std::valarray<T>* v;
	std::slice s;
	size_type curr;	// index of current element

	T& ref(size_type i) const { return (*v)[s.start()+i*s.stride()]; }
public:
	// no copy in constructors
	Slice_iter(std::valarray<T>* vv, std::slice ss) :v(vv), s(ss), curr(0) { }
	Slice_iter(const Slice_iter& it) : v(it.v), s(it.s), curr(it.curr) {}
	// copy in assignment
	Slice_iter& operator=(const Slice_iter& it);

	Slice_iter begin() const;
	Slice_iter end() const;
	size_type size() const { return s.size(); }

	Slice_iter& operator++() { ++curr; return *this; }
	Slice_iter operator++(int) { Slice_iter it = *this; ++curr; return it; }

	T& operator[](size_type i) const { return ref(i); }		// C style subscript
	T& operator()(size_type i) const { return ref(i); }		// Fortran-style subscript
	T& operator*() { return ref(curr); }			// current element
	T* operator->() { T& curr_elem = ref(curr); return &curr_elem; }

	friend bool operator==<>(const Slice_iter& it_a, const Slice_iter& it_b);
	friend bool operator!=<>(const Slice_iter& it_a, const Slice_iter& it_b);
	friend bool operator< <>(const Slice_iter& it_a, const Slice_iter& it_b);

	template<typename T2>
	friend std::ostream& operator<<(std::ostream& os, const Slice_iter<T2>& q);
};


template<typename T>
Slice_iter<T>& Slice_iter<T>::operator=(const Slice_iter<T>& it) 
{ 
	if(&it != this)
	{
		assert(s.size() == it.s.size());
		for(size_type i = 0; i<s.size(); ++i)
			ref(i) = it.ref(i);
		s = it.s;
		curr = it.curr;
	} 
	return *this;
}


template<typename T>
Slice_iter<T> Slice_iter<T>::begin() const
{
	return Slice_iter(v, s);
}


template<typename T>
Slice_iter<T> Slice_iter<T>::end() const
{
	Slice_iter t = *this;
	t.curr = s.size();	// index of last-plus-one element
	return t;
}


template<typename T>
bool operator==(const Slice_iter<T>& it_a, const Slice_iter<T>& it_b)
{
	return it_a.curr==it_b.curr
		&& it_a.s.stride()==it_b.s.stride()
		&& it_a.s.start()==it_b.s.start();
}


template<typename T>
bool operator!=(const Slice_iter<T>& it_a, const Slice_iter<T>& it_b)
{
	return !(it_a==it_b);
}

template<typename T>
bool operator<(const Slice_iter<T>& it_a, const Slice_iter<T>& it_b)
{
	return it_a.curr<it_b.curr
		&& it_a.s.stride()==it_b.s.stride()
		&& it_a.s.start()==it_b.s.start();
}


template<typename T>
std::ostream& operator<<(std::ostream& os, const Slice_iter<T>& it)
{
	int i=0;
	os << '(';
	for(Slice_iter<T> it2(it); it2 != it2.end(); ++it2, ++i)
	{
		std::cout << *it2;
		if(i != it.size()-1) os << " ";
	} 
	os << ')';
	return os;
}


template<typename T>
std::valarray<T> scale_and_add(const std::valarray<T>& a,
				const T& b, const std::valarray<T>& c)
{
	return ((a * b) + c);
}


template<typename T>
Slice_iter<T> scale_and_add(const Slice_iter<T>& it_a, const T& b, const Slice_iter<T>& it_c)
{
	Slice_iter<T> it_res(it_c);
	for(Slice_iter<T> it_a2(it_a); it_a2 != it_a2.end(); ++it_res, ++it_a2)
		*it_res += *it_a2 * b;
	return it_res;
}


template<typename T>
T dot_product(const Slice_iter<T>& it_a, const Slice_iter<T>& it_b)
{
	T res = 0.;
	for(Slice_iter<T> it_a2(it_a), it_b2(it_b); it_a2 != it_a2.end(); ++it_a2, ++it_b2)
		res += *it_a2 * *it_b2;

}

//-------------------------------------------------------------



// forward declarations to allow friend declarations:
template<typename T> class Cslice_iter;
template<typename T> bool operator==(const Cslice_iter<T>&, const Cslice_iter<T>&);
template<typename T> bool operator!=(const Cslice_iter<T>&, const Cslice_iter<T>&);
template<typename T> bool operator< (const Cslice_iter<T>&, const Cslice_iter<T>&);


template<typename T> class Cslice_iter : public std::iterator<std::forward_iterator_tag, T>
{
	const std::valarray<T>* v;
	std::slice s;
	size_type curr; // index of current element
	const T& ref(size_type i) const { return (*v)[s.start()+i*s.stride()]; }
public:
	Cslice_iter(const std::valarray<T>* vv, std::slice ss): v(vv), s(ss), curr(0){}
	Cslice_iter(const Cslice_iter& cit) : v(cit.v),  s(cit.s), curr(cit.curr) {}
	Cslice_iter end() const
	{
		Cslice_iter cit = *this;
		cit.curr = s.size(); // index of one plus last element
		return cit;
	}
	Cslice_iter& operator++() { curr++; return *this; }
	Cslice_iter operator++(int) { Cslice_iter cit = *this; ++curr; return cit; }
	
	const T& operator[](size_type i) const { return ref(i); }
	const T& operator()(size_type i) const { return ref(i); }
	const T& operator*() const { return ref(curr); }
	const T* operator->() const { const T& curr_elem = ref(curr); return &curr_elem; }

	friend bool operator==<>(const Cslice_iter& cit_a, const Cslice_iter& cit_b);
	friend bool operator!=<>(const Cslice_iter& cit_a, const Cslice_iter& cit_b);
	friend bool operator< <>(const Cslice_iter& cit_a, const Cslice_iter& cit_b);

	template<typename T2>
	friend std::ostream& operator<<(std::ostream& os, const Cslice_iter<T2>& cit);
};


template<typename T>
bool operator==(const Cslice_iter<T>& cit_a, const Cslice_iter<T>& cit_b)
{
	return cit_a.curr==cit_b.curr
		&& cit_a.s.stride()==cit_b.s.stride()
		&& cit_a.s.start()==cit_b.s.start();
}


template<typename T>
bool operator!=(const Cslice_iter<T>& cit_a, const Cslice_iter<T>& cit_b)
{
	return !(cit_a==cit_b);
}


template<typename T>
bool operator<(const Cslice_iter<T>& cit_a, const Cslice_iter<T>& cit_b)
{
	return cit_a.curr<cit_b.curr
		&& cit_a.s.stride()==cit_b.s.stride()
		&& cit_a.s.start()==cit_b.s.start();
}


template<typename T>
std::ostream& operator<<(std::ostream& os, const Cslice_iter<T>& cit)
{
	int i=0;
	os << '(';
	for(Cslice_iter<T> cit2(cit); cit2<cit2.end(); ++cit2, ++i) // it < q.end()
	{
		std::cout << *cit2;
		if(i != cit.s.size()-1) os << " ";
	} 
	os << ')';
	return os;
}


template<typename T>
T dot_product(const Cslice_iter<T>& cit_a, const Cslice_iter<T>& cit_b)
{
	T res = 0.;
	for(Cslice_iter<T> cit_a2(cit_a), cit_b2(cit_b); cit_a2 != cit_a.end(); ++cit_a2, ++cit_b2)
		res += *cit_a2 * *cit_b2;
	return res;
}

#endif
