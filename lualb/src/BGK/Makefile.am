SUBDIRS = . 

#
#  main sources
#

lbBGKsources = 			\
	lb-BGK.h		\
	lb-BGK.c		\
	lb-d2q9-BGK.h		\
	lb-d2q9-BGK.c		\
	lb-d3q19-BGK.h		\
	lb-d3q19-BGK.c

noinst_LIBRARIES = liblbBGK.a

liblbBGK_a_SOURCES = 		\
	$(lbBGKsources)	
	
liblbBGK_a_LIBDADD = ../core/liblbcore.a $(HDF5_LIBS)

INCLUDES = -I../core -I$(top_srcdir) $(HDF5_CFLAGS)

EXTRA_DIST = README
