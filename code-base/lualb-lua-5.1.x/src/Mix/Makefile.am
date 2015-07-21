SUBDIRS = . 

#
#  main sources
#

lbMixsources = 			\
	lb-Mix.h		\
	lb-Mix.c		\
	lb-d2q9-Mix.h		\
	lb-d2q9-Mix.c		\
	lb-d3q19-Mix.h		\
	lb-d3q19-Mix.c

noinst_LIBRARIES = liblbMix.a

liblbMix_a_SOURCES = 		\
	$(lbMixsources)	
	
liblbMix_a_LIBDADD = ../core/liblbcore.a $(HDF5_LIBS)

INCLUDES = -I../core -I$(top_srcdir) $(HDF5_CFLAGS)

EXTRA_DIST = README
