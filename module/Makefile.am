##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(srcdir) \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

lib_LTLIBRARIES = \
  libcpm.la

libcpm_la_LIBADD = \
  -lCore \
  -lffaobjects \
  -lphool \
  -lRIO \
  -lSubsysReco \
  -lTree \
  -ltrack_io \
  -ltrackbase_historic_io \
  -ltpc_io

pkginclude_HEADERS = \
  CPMLocalLinePoCA.h \
  CPMRecord.h \
  CPMTypes.h \
  CPMVoxelContainer.h \
  PHCPMTpcCalibration.h

libcpm_la_SOURCES = \
  CPMLocalLinePoCA.cc \
  PHCPMTpcCalibration.cc

noinst_PROGRAMS = \
  test_local_line_poca

test_local_line_poca_SOURCES = \
  test_local_line_poca.cc

test_local_line_poca_LDADD = \
  libcpm.la

##############################################
