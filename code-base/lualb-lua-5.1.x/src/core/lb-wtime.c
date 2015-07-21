#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <stdio.h>

#ifdef LB_ENABLE_MPI
# include <mpi.h>
#else /* !LB_ENABLE_MPI */
# include <sys/time.h>
#endif /* LB_ENABLE_MPI */

#include "lb-wtime.h"

double
lb_wtime(void)
{
#ifdef LB_ENABLE_MPI
	return MPI_Wtime();
#else /* !LB_ENABLE_MPI */
	struct timeval tv;
	struct timezone tz;

	gettimeofday(&tv, &tz);

	return tv.tv_sec + tv.tv_usec/1.0e6;
#endif /* LB_ENABLE_MPI */
}

const char*
lb_wtime_string(double sec)
{
#define BUF_SIZE 64
	static char buf[BUF_SIZE + 1];

	if (sec/60 < 1.0)
		snprintf(buf, BUF_SIZE, "%f sec", sec);
	else if (sec/(60*60) < 1.0) {
		const double min = sec/60;
		snprintf(buf, BUF_SIZE, "%f min", min);
	} else {
		const double hours = sec/(60*60);
		snprintf(buf, BUF_SIZE, "%f hours", hours);
	}

	return buf;
}
