#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef LB_ENABLE_MPI
#include <mpi.h>
#endif /* LB_ENABLE_MPI */

#include "lb-messages.h"

static int lb_info_is_enabled = TRUE;

#ifdef LB_ENABLE_MPI
static int
mpi_initialized(void)
{
	int flag;

	MPI_Initialized(&flag);

	return flag;
}
#endif /* LB_ENABLE_MPI */

void
lb_abort(void)
{
#ifdef LB_ENABLE_MPI
	if (mpi_initialized())
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif /* LB_ENABLE_MPI */
	abort();
}

void
lb_fatal(const char* fmt, ...)
{
	va_list ap;

#ifdef LB_ENABLE_MPI
	int processor_name_length;
	char processor_name[MPI_MAX_PROCESSOR_NAME] = "<MPI not initialized>";

	if (mpi_initialized())
		MPI_Get_processor_name(processor_name, &processor_name_length);

	fprintf(stderr, "** ERROR ** (%s:%d): ",
		processor_name, getpid());
#else
	fprintf(stderr, "** ERROR ** (process:%d): ", getpid());
#endif /* LB_ENABLE_MPI */

	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);

	fflush(stderr);

	lb_abort();
}

void
lb_warning(const char* fmt, ...)
{
	va_list ap;

	lb_assert(fmt != NULL);

#ifdef LB_ENABLE_MPI
	int processor_name_length;
	char processor_name[MPI_MAX_PROCESSOR_NAME] = "<MPI not initialized>";

	if (mpi_initialized())
		MPI_Get_processor_name(processor_name, &processor_name_length);

	fprintf(stderr, "** WARNING ** (%s:%d): ",
		processor_name, getpid());
#else
	fprintf(stderr, "** WARNING ** (process:%d): ", getpid());
#endif /* LB_ENABLE_MPI */

	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);

	fflush(stderr);
}

void
lb_info(const char* fmt, ...)
{
	va_list ap;

	lb_assert(fmt != NULL);

	if (lb_info_is_enabled) {
		va_start(ap, fmt);
		vfprintf(stdout, fmt, ap);
		va_end(ap);
		fflush(stdout);
	}
}

void
lb_info_enable(void)
{
	lb_info_is_enabled = TRUE;
}

void
lb_info_disable(void)
{
	lb_info_is_enabled = FALSE;
}

int
lb_info_enabled(void)
{
	return lb_info_is_enabled;
}

int
lb_is_parallel(void)
{

#ifdef LB_ENABLE_MPI
return 1;
#endif /* LB_ENABLE_MPI */

return 0;
}
