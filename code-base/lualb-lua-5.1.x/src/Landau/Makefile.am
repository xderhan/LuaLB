SUBDIRS = . 

#
#  main sources
#

lbLandausources = 			\
	lb-Landau.h		\
	lb-Landau.c		\
	lb-d2q9-Landau.h		\
	lb-d2q9-Landau.c		\
	lb-d3q19-Landau.h		\
	lb-d3q19-Landau.c

noinst_LIBRARIES = liblbLandau.a

liblbLandau_a_SOURCES = 		\
	$(lbLandausources)	
	
liblbLandau_a_LIBDADD = ../core/liblbcore.a $(HDF5_LIBS)

INCLUDES = -I../core -I$(top_srcdir) $(HDF5_CFLAGS)

EXTRA_DIST = README
