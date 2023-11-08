#include <errno.h>
#include <filesystem>
#include <iostream>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/statvfs.h>

#include "LowTasks.h"
#include "backtrace.h"
#include "global.h"
#include "logging.h"

static void _mkdir(const char *dir, mode_t mode)
{
    // from
    // https://stackoverflow.com/questions/2336242/recursive-mkdir-system-call-on-unix
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf(tmp, sizeof(tmp), "%s", dir);
    len = strlen(tmp);
    if (tmp[len - 1] == '/')
	tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++)
	if (*p == '/') {
	    *p = 0;
	    mkdir(tmp, mode);
	    *p = '/';
	}
    const int res = mkdir(tmp, S_IRWXU);
    if (res != 0) {
	logging::print_master(LOG_ERROR
			      "mkdir returned %d for path %s and mode %d\n",
			      res, tmp, mode);
    }
}

void ensure_directory_exists(const std::string &dirname)
{
    if (CPU_Master) {
	struct stat buffer;
	if (stat(dirname.c_str(), &buffer)) {
	    _mkdir(dirname.c_str(), 0755);
	}
    }
}

void delete_directory_if_exists(const std::string &dirname)
{
    if (CPU_Master) {
	if (std::filesystem::exists(dirname)) {
	    std::filesystem::remove_all(dirname);
	}
    }
}

/**
	Finalize MPI and terminate program.

	\param returncode returncode to return
*/
void PersonalExit(int returncode)
{
    std::flush(std::cout);
    if (returncode != 0) {
	PrintTrace();
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
    fprintf(stderr, "\n\nfatal: %s\n", msg);
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

std::string lowercase(const std::string& s)
{
    std::string out;
    out.resize(s.size());
    std::transform(s.begin(), s.end(),
		   out.begin(), ::tolower);
	return out;
}
