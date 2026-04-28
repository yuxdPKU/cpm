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
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(OPT_SPHENIX)/lib

lib_LTLIBRARIES = \
  libcpm.la

libcpm_la_LIBADD = \
  -lActsCore \
  -lCore \
  -lffaobjects \
  -lphool \
  -lRIO \
  -lSubsysReco \
  -lTree \
  -ltpc \
  -ltrack_io \
  -ltrack \
  -ltrackbase_historic_io \
  -ltrackbase_historic \
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
  testexternals_cpm \
  test_local_line_poca

BUILT_SOURCES = \
  testexternals.cc

testexternals_cpm_SOURCES = \
  testexternals.cc

testexternals_cpm_LDADD = \
  libcpm.la

test_local_line_poca_SOURCES = \
  CPMLocalLinePoCA.cc \
  test_local_line_poca.cc

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)

##############################################
