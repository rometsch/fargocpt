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
#include "parameters.h"

#include "output.h"
void handle_sigterm_outputs(t_data &data)
{
    logging::print_master(
	LOG_INFO "Received SIGTERM, starting to writing debug output\n");
    // Enable output of Qplus / Qminus for bitwise exact restarting.
    if (!data[t_data::QPLUS].get_write()) {
	data[t_data::QPLUS].set_write(true, false);
    }
    if (!data[t_data::QMINUS].get_write()) {
	data[t_data::QMINUS].set_write(true, false);
    }

    if (parameters::variableGamma) {
	if (!data[t_data::GAMMAEFF].get_write()) {
	    data[t_data::GAMMAEFF].set_write(true, false);
	    data[t_data::MU].set_write(true, false);
	    data[t_data::GAMMA1].set_write(true, false);
	}
    }

    output::write_grids(data, N_output, N_hydro_iter, PhysicalTime, true);
    data.get_planetary_system().write_planets(N_output, 2);
    output::write_misc(true);

    MPI_Barrier(MPI_COMM_WORLD);
    logging::print_master(LOG_INFO "Wrote debug outputs\n");
    die("Received SIGTERM\n");
}

/**
	Finalize MPI and terminate program.

	\param returncode returncode to return
*/
void PersonalExit(int returncode)
{
    std::flush(std::cout);
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
    PersonalExit(128);
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
