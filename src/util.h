#ifndef _UTIL_H_
#define _UTIL_H_

template< typename T >
inline T pow2(T x)
{
	return x*x;
}

template< typename T >
inline T pow3(T x)
{
	return x*x*x;
}

template< typename T >
inline T pow4(T x)
{
	return x*x*x*x;
}

template< typename T >
inline T min(T x, T y)
{
	return x < y ? x : y;
}

template< typename T >
inline T max(T x, T y)
{
	return x > y ? x : y;
}

double cutoff(double point, double width, double x);

#endif
