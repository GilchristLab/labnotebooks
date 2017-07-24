cmake_minimum_required(VERSION 3.5)
project(src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D STANDALONE")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D DENIZHAN")

set(SOURCE_FILES
    include/base/Model.h
    include/base/Parameter.h
    include/base/Trace.h
    include/FONSE/FONSEModel.h
    include/FONSE/FONSEParameter.h
    include/PA/PAModel.h
    include/PA/PAParameter.h
    include/ROC/ROCModel.h
    include/ROC/ROCParameter.h
    include/CodonTable.h
    include/CovarianceMatrix.h
    include/Gene.h
    include/Genome.h
    include/MCMCAlgorithm.h
    include/mixtureDefinition.h
    include/SequenceSummary.h
    include/Testing.h
    include/Utility.h
    CodonTable.cpp
    CovarianceMatrix.cpp
    FONSEModel.cpp
    FONSEParameter.cpp
    Gene.cpp
    Genome.cpp
    main.cpp
    MCMCAlgorithm.cpp
    Model.cpp
    Parameter.cpp
    RCPP_Model.cpp
    RCPP_Parameter.cpp
    RCPP_Trace.cpp
    PAModel.cpp
    PAParameter.cpp
    ROCModel.cpp
    ROCParameter.cpp
    SequenceSummary.cpp
    Testing.cpp
    Trace.cpp)

add_executable(src ${SOURCE_FILES})
