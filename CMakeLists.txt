cmake_minimum_required (VERSION 2.8)

project(BAMStats CXX)

  # set our library and executable destination dirs
  set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )

  add_executable(BAMStats Bam_Stats.cpp alglibinternal.cpp alglibmisc.cpp ap.cpp dataanalysis.cpp diffequations.cpp fasttransforms.cpp integration.cpp interpolation.cpp linalg.cpp optimization.cpp solvers.cpp specialfunctions.cpp statistics.cpp continued_fraction.cpp load_data_for_complexity.cpp moment_sequence.cpp OptionParser.cpp SAM.cpp RNG.cpp QualityScore.cpp MappedRead.cpp GenomicRegion.cpp FileIterator.cpp smithlab_utils.cpp smithlab_os.cpp sim_utils.cpp)

  set_property(TARGET BAMStats PROPERTY CXX_STANDARD 11)
  set_property(TARGET BAMStats PROPERTY CXX_STANDARD_REQUIRED ON)
 
  find_package(ZLIB REQUIRED)
  add_library(libbamtools STATIC IMPORTED)
  add_library(libbam STATIC IMPORTED)
  add_library(libgsl STATIC IMPORTED)
  add_library(libgslcblas STATIC IMPORTED)
  set_target_properties(libbamtools PROPERTIES IMPORTED_LOCATION /sw/apps/bioinfo/bamtools/2.3.0/milou/lib/libbamtools.so.2.3.0)        
  set_target_properties(libbam PROPERTIES IMPORTED_LOCATION /sw/apps/bioinfo/samtools/0.1.9/libbam.a)
  set_target_properties(libgsl PROPERTIES IMPORTED_LOCATION /sw/apps/gsl/1.16/milou/lib/libgsl.so)        
  set_target_properties(libgslcblas PROPERTIES IMPORTED_LOCATION /sw/apps/gsl/1.16/milou/lib/libgslcblas.so)        
  include_directories(/pica/h1/pelin/External/boost_1_58_0)
  include_directories(/sw/apps/bioinfo/samtools/0.1.9)
  include_directories(/pica/h1/pelin/bin/Bam_stats/preseq-master)
  include_directories(/pica/h1/pelin/bin/Bam_stats/preseq-master/smithlab_cpp)
  include_directories(/sw/apps/bioinfo/bamtools/2.3.0/milou/include)
  link_directories(/sw/apps/gsl/1.16/milou/lib)
  link_directories(/sw/apps/bioinfo/samtools/0.1.9)
  link_directories(/sw/apps/bioinfo/bamtools/2.3.0/milou/lib)
  target_link_libraries(BAMStats libbamtools libgsl libgslcblas libbam ${ZLIB_LIBRARIES})
