/**
	\file radialarray.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>

	This file manages the t_radialarray class which stores radial
   information
*/

#include "radialarray.h"
#include <stdlib.h>

t_radialarray::t_radialarray()
{
    m_size = 0;
    array = NULL;
}

t_radialarray::t_radialarray(ptrdiff_t size)
{
    // assign memory
    m_size = size;
    array = new double[m_size];

    // set to 0
    clear();
}

t_radialarray::~t_radialarray() { delete[] array; }

/**
	set all entries to 0
*/
void t_radialarray::clear()
{
    for (ptrdiff_t i = 0; i < m_size; ++i) {
	operator()(i) = 0;
    }
}

/**
	resize radial array. radial is set to 0 automatically
*/
void t_radialarray::resize(ptrdiff_t size)
{
    // delete old data
    delete[] array;

    // assign new memory
    m_size = size;
    array = new double[m_size];

    // set to 0
    clear();
}
