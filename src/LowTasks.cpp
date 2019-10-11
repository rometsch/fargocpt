#include <errno.h>
#include <iostream>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "LowTasks.h"
#include "global.h"
#include "logging.h"

double GetGlobalIFrac(double r)
{
	int i = 0;
	double ifrac;

	if (r < GlobalRmed[0])
		return 0.0;
	if (r > GlobalRmed[GlobalNRadial - 1])
		return (double)GlobalNRadial - 1.0;
	while (GlobalRmed[i] <= r)
		i++;
	ifrac = (double)i +
		(r - GlobalRmed[i - 1]) / (GlobalRmed[i] - GlobalRmed[i - 1]) -
		1.0;

	return ifrac;
}

/**
	Finalize MPI and terminate program.

	\param returncode returncode to return
*/
void PersonalExit(int returncode)
{
	std::flush(std::cout);
	MPI_Barrier(MPI_COMM_WORLD);
	if (returncode != 0) {
		MPI_Abort(MPI_COMM_WORLD, returncode);
	}
	MPI_Finalize();
	exit(returncode);
}

t_polargrid *CreatePolarGrid(unsigned int Nr, unsigned int Ns, const char *name)
{
	t_polargrid *polarGrid;

	polarGrid = new t_polargrid();
	polarGrid->set_size(Nr, Ns);
	polarGrid->set_name(name);
	polarGrid->WriteOut = 0;
	polarGrid->Initialised = 0;

	return polarGrid;
}

void MultiplyPolarGridbyConstant(t_polargrid *arraysrc, double constant)
{
	(*arraysrc) *= constant;
}

static void die_builtin(const char *err, va_list params)
{
	char msg[1024];
	vsnprintf(msg, sizeof(msg), err, params);
	fprintf(stderr, "fatal: %s\n", msg);
	MPI_Finalize();
	exit(128);
}

void die(const char *err, ...)
{
	va_list params;

	va_start(params, err);
	die_builtin(err, params);
	va_end(params);
}

void die_errno(const char *fmt, ...)
{
	va_list params;
	char fmt_with_err[1024];
	char str_error[256], *err;
	unsigned int i, j;

	err = strerror(errno);
	for (i = j = 0; err[i] && j < sizeof(str_error) - 1;) {
		if ((str_error[j++] = err[i++]) != '%')
			continue;
		if (j < sizeof(str_error) - 1) {
			str_error[j++] = '%';
		} else {
			/* No room to double the '%', so we overwrite it with
			 * '\0' below */
			j--;
			break;
		}
	}
	str_error[j] = 0;
	snprintf(fmt_with_err, sizeof(fmt_with_err), "%s: %s", fmt, str_error);

	va_start(params, fmt);
	die_builtin(fmt_with_err, params);
	va_end(params);
}
