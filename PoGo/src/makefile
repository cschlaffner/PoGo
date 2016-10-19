CXX = g++
CXXFLAGS = -Wall -g -std=c++0x 


all: main clean

main: main.o EnumStringMapper.o Utils.o FastaEntry.o FastaParser.o ProteinEntry.o CoordinateWrapper.o PeptideCoordinates.o GeneEntry.o PTMEntry.o PeptideEntry.o KmerMap.o PeptideMapperUtils.o PossibleKeyGenerator.o MapEntry.o MappedPeptides.o ResultParser.o GTFParser.o Globals.o
	$(CXX) $(CXXFLAGS) -o PoGo main.o EnumStringMapper.o Utils.o FastaEntry.o FastaParser.o ProteinEntry.o CoordinateWrapper.o PeptideCoordinates.o GeneEntry.o PTMEntry.o PeptideEntry.o KmerMap.o PossibleKeyGenerator.o MapEntry.o MappedPeptides.o ResultParser.o GTFParser.o PeptideMapperUtils.o Globals.o

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp 

GTFParser.o: GTFParser.h ResultParser.h

ResultParser.o: ResultParser.h MappedPeptides.h

MappedPeptides.o: MappedPeptides.h MapEntry.h

MapEntry.o: MapEntry.h KmerMap.h

KmerMap.o: KmerMap.h KmerEntry.h PeptideEntry.h PeptideMapperUtils.h PossibleKeyGenerator.h

PeptideEntry.o: PeptideEntry.h PTMEntry.h

PeptideMapperUtils.o: PeptideMapperUtils.h Utils.h

PossibleKeyGenerator.o: PossibleKeyGenerator.h PeptideMapperUtils.h

PTMEntry.o: PTMEntry.h GeneEntry.h

GeneEntry.o: GeneEntry.h PeptideCoordinates.h

PeptideCoordinates.o: PeptideCoordinates.h CoordinateWrapper.h

CoordinateWrapper.o: CoordinateWrapper.h ProteinEntry.h

ProteinEntry.o: ProteinEntry.h FastaParser.h

FastaParser.o: FastaParser.h FastaEntry.h 

FastaEntry.o: FastaEntry.h Utils.h

Utils.o: Utils.h EnumStringMapper.h

EnumStringMapper.o: EnumStringMapper.h Coordinates.h Globals.h

Globals.o: Globals.h

clean:
	rm -f *.o

