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

ACTS_LIBS = \
  -lActsCore \
  -lActsPluginTGeo \
  -lActsExamplesDetectorTGeo \
  -lActsExamplesFramework

CPM_EXTERNAL_LIBS = \
  $(ACTS_LIBS) \
  -lCore \
  -lffaobjects \
  -lHist \
  -lphool \
  -lPhysics \
  -lRIO \
  -lSubsysReco \
  -lTree \
  -ltpc \
  -ltrack_io \
  -ltrack \
  -ltrackbase_historic_io \
  -ltrackbase_historic \
  -ltpc_io \
  -ltpccalib

libcpm_la_LIBADD = \
  $(CPM_EXTERNAL_LIBS)

pkginclude_HEADERS = \
  CPMAverageCorrectionReconstruction.h \
  CPMRecord.h \
  CPMReconstructionHelper.h \
  CPMVoxelContainer.h \
  CPMVoxelContainerv1.h \
  PHCPMTpcCalibration.h \
  PHCPMTpcCalibrationQA.h

ROOTDICTS = \
  CPMVoxelContainer_Dict.cc \
  CPMVoxelContainerv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  CPMVoxelContainer_Dict_rdict.pcm \
  CPMVoxelContainerv1_Dict_rdict.pcm

libcpm_la_SOURCES = \
  $(ROOTDICTS) \
  CPMAverageCorrectionReconstruction.cc \
  CPMReconstructionHelper.cc \
  CPMVoxelContainerv1.cc \
  PHCPMTpcCalibration.cc \
  PHCPMTpcCalibrationQA.cc

noinst_PROGRAMS = \
  testexternals_cpm \
  test_local_line_poca

BUILT_SOURCES = \
  testexternals.cc

testexternals_cpm_SOURCES = \
  testexternals.cc

testexternals_cpm_LDADD = \
  libcpm.la \
  $(CPM_EXTERNAL_LIBS)

test_local_line_poca_SOURCES = \
  test_local_line_poca.cc

test_local_line_poca_LDADD = \
  libcpm.la \
  $(CPM_EXTERNAL_LIBS)

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

# just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)

##############################################
