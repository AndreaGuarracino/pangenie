add_library(PanGenieLib SHARED 
	emissionprobabilitycomputer.cpp
	copynumber.cpp
	commandlineparser.cpp
	commands.cpp
	columnindexer.cpp
	dnasequence.cpp
	fastareader.cpp
	genotypingresult.cpp
	graphbuilder.cpp
	graph.cpp
	haplotypesampler.cpp
	histogram.cpp
	hmm.cpp
	jellyfishcounter.cpp
	staticdictcounter.cpp
	sshash_backend.cpp
	jellyfishreader.cpp
	kmerpath.cpp
	kmerpath16.cpp
	kmerparser.cpp
	pathsampler.cpp
	probabilitycomputer.cpp
	probabilitytable.cpp
	sampledpanel.cpp
	samplingemissions.cpp
	samplingtransitions.cpp
	sequenceutils.cpp
	stepwiseuniquekmercomputer.cpp
	timer.cpp
	transitionprobabilitycomputer.cpp
	threadpool.cpp
	uniquekmercomputer.cpp
	multiallelicuniquekmers.cpp
	biallelicuniquekmers.cpp
	variant.cpp
	variantreader.cpp)

# Edit D1: SSHash static-dictionary backend. Only sshash_backend.cpp needs these headers.
# SSHASH_ROOT points at an sshash source checkout (github.com/jermp/sshash, recursive).
# Default is a repo-relative external/sshash; override with -DSSHASH_ROOT=/path.
if(NOT DEFINED SSHASH_ROOT)
	set(SSHASH_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../external/sshash")
endif()
target_include_directories(PanGenieLib PRIVATE
	${SSHASH_ROOT}
	${SSHASH_ROOT}/include
	${SSHASH_ROOT}/external/pthash/include
	${SSHASH_ROOT}/external/pthash/external/bits/include
	${SSHASH_ROOT}/external/pthash/external/fastmod
	${SSHASH_ROOT}/external/pthash/external/mm_file/include
	${SSHASH_ROOT}/external/pthash/external/cmd_line_parser/include
	${SSHASH_ROOT}/external/pthash/external/xxHash
	${SSHASH_ROOT}/external/pthash/external/bits/external/essentials/include
	${SSHASH_ROOT}/external/cityhash
	${SSHASH_ROOT}/external/gz)


add_executable(PanGenie-index pangenie-index.cpp)
target_link_libraries(PanGenie-index PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-index PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})


add_executable(PanGenie pangenie-genotype.cpp)
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})


add_executable(PanGenie-vcf pangenie-vcf.cpp)
target_link_libraries(PanGenie-vcf PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-vcf PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})


add_executable(PanGenie-sampling pangenie-sampling.cpp)
target_link_libraries(PanGenie-sampling PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(PanGenie-sampling PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})


add_executable(Analyze-UK analyze-uk.cpp)
target_link_libraries(Analyze-UK PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(Analyze-UK PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})
