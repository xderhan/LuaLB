SUBDIRS = . 

#
#  graphics sources
#

graphicssources = 

lbrgbsources =			\
	lb-rgb.h		\
	lb-rgb.c		\
	lb-colormap.h		\
	lb-colormap.c

if LB_ENABLE_RGB
maybe_lbrgbsources = $(lbrgbsources) 
else
maybe_lbrgbsources =
endif

#
# library
#

noinst_LIBRARIES = liblbgraphics.a

liblbgraphics_a_SOURCES = 	\
	$(lbutilssources)	\
	$(maybe_lbrgbsources)

liblbgraphics_a_LIBDADD = $(PNG_LIBS)

EXTRA_liblbgraphics_a_SOURCES = \
	$(lbrgbsources)		

INCLUDES = -I../core -I$(top_srcdir)

EXTRA_DIST = README 
