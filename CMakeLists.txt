cmake_minimum_required(VERSION 3.6)
project(PoGo)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    PoGo/src/Coordinates.h
    PoGo/src/CoordinateWrapper.cpp
    PoGo/src/CoordinateWrapper.h
    PoGo/src/EnumStringMapper.cpp
    PoGo/src/EnumStringMapper.h
    PoGo/src/ExistingPeptides.h
    PoGo/src/FastaEntry.cpp
    PoGo/src/FastaEntry.h
    PoGo/src/FastaParser.cpp
    PoGo/src/FastaParser.h
    PoGo/src/GeneEntry.cpp
    PoGo/src/GeneEntry.h
    PoGo/src/Globals.cpp
    PoGo/src/Globals.h
    PoGo/src/GTFParser.cpp
    PoGo/src/GTFParser.h
    PoGo/src/KmerEntry.h
    PoGo/src/KmerMap.cpp
    PoGo/src/KmerMap.h
    PoGo/src/main.cpp
    PoGo/src/main_functions.h
    PoGo/src/MapEntry.cpp
    PoGo/src/MapEntry.h
    PoGo/src/MappedPeptides.cpp
    PoGo/src/MappedPeptides.h
    PoGo/src/PeptideCoordinates.cpp
    PoGo/src/PeptideCoordinates.h
    PoGo/src/PeptideEntry.cpp
    PoGo/src/PeptideEntry.h
    PoGo/src/PeptideMapperUtils.cpp
    PoGo/src/PeptideMapperUtils.h
    PoGo/src/PossibleKeyGenerator.cpp
    PoGo/src/PossibleKeyGenerator.h
    PoGo/src/ProteinEntry.cpp
    PoGo/src/ProteinEntry.h
    PoGo/src/PTMEntry.cpp
    PoGo/src/PTMEntry.h
    PoGo/src/ResultParser.cpp
    PoGo/src/ResultParser.h
    PoGo/src/Utils.cpp
    PoGo/src/Utils.h
    LICENSE
    README.md)

add_executable(PoGo ${SOURCE_FILES})