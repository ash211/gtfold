INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

AM_CFLAGS = $(OPENMP_CFLAGS) -DDATADIR='$(datadir)/@PACKAGE@'
AM_CXXFLAGS = -DDATADIR='$(datadir)/@PACKAGE@'

bin_PROGRAMS = gtfold

gtfold_SOURCES = \
	main.cc\
	loader.cc\
	algorithms.c\
	algorithms-partition.c\
	traceback.c


gtfold_LDFLAGS = 

gtfold_LDADD = -lm

CLEANFILES = *~ *.o
