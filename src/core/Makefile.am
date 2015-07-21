SUBDIRS = . 

#
#  core sources
#

lbsources = 			\
	lb-macros.h		\
	lb-partition-info.h	\
	lb-memory.h		\
	lb-memory.c		\
	lb-messages.h		\
	lb-messages.c		\
	lb-wtime.h		\
	lb-wtime.c		\
	lb-d2q9.h		\
	lb-d2q9.c		\
	lb-d3q19.h		\
	lb-d3q19.c		\
	lb-d2q9-utils.h  	\
	lb-d2q9-utils.c  	\
	lb-d3q19-utils.h 	\
	lb-d3q19-utils.c 

lbmpisources =			\
	lb-mpi-utils.h		\
	lb-mpi-utils.c		\
	lb-cart-partition.h	\
	lb-cart-partition.c

if LB_ENABLE_MPI
maybe_lbmpisources = $(lbmpisources)
else
maybe_lbmpisources =
endif

#
# library
#

noinst_LIBRARIES = liblbcore.a

liblbcore_a_SOURCES = 		\
	$(lbsources)		\
	$(maybe_lbmpisources)	

EXTRA_liblbcore_a_SOURCES = 	\
	$(lbmpisources)		

INCLUDES = -I$(top_srcdir) 

EXTRA_DIST = README 
