set (CMAKE_CXX_STANDARD 20)
set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
include_directories (${PROGRAM_SOURCE_DIR})
add_executable(tests tests.cpp utils.cpp EmissionProbabilityComputerTest.cpp CopyNumberTest.cpp UniqueKmersTest.cpp UniqueKmerComputerTest.cpp KmerPathTest.cpp VariantTest.cpp VariantReaderTest.cpp GraphBuilderTest.cpp ProbabilityComputerTest.cpp TransitionProbabilityComputerTest.cpp HMMTest.cpp ColumnIndexerTest.cpp GenotypingResultTest.cpp DnaSequenceTest.cpp FastaReaderTest.cpp KmerCounterTest.cpp HistogramTest.cpp PathSamplerTest.cpp ProbabilityTableTest.cpp KmerParser.cpp HaplotypeSamplerTest.cpp SamplingEmissionsTest.cpp SamplingTransitionsTest.cpp SampledPanelTest.cpp CommandsTest.cpp)

target_link_libraries(tests PanGenieLib ${JELLYFISH_LDFLAGS_OTHER} ${ZLIB_LDFLAGS_OTHER} ${CEREAL_LDFLAGS_OTHER})
target_link_libraries(tests PanGenieLib ${JELLYFISH_LIBRARIES} ${ZLIB_LIBRARIES} ${CEREAL_LIBRARIES})

# Tests use paths relative to their working directory and create output
# fixtures. Keep those writes out of the source tree.
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/data" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
add_test(NAME pangenie-tests COMMAND tests)
set_tests_properties(pangenie-tests PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
