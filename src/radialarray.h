#pragma once

#include <stddef.h>

class t_radialarray
{
  private:
    /// number of entries (size of array)
    ptrdiff_t m_size;

  public:
    t_radialarray();
    t_radialarray(ptrdiff_t size);
    ~t_radialarray();

    /// pointer to the array
    double *array;

    void clear();
    void resize(ptrdiff_t size);

    inline double &operator()(ptrdiff_t n_radial) { return array[n_radial]; }
    inline double operator()(ptrdiff_t n_radial) const
    {
	return array[n_radial];
    }

    inline operator const double *() const { return array; }
    inline operator double *() { return array; }
};
