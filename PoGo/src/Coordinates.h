#ifndef COORDINATES_H
#define COORDINATES_H

#include "Globals.h"

//possible offsets.
enum Offset {
	off1 = 1,
	off2 = 2,
	off3 = 3
};

//possible chromosomes
enum Chromosome {
	chr1 = 1,
	chr1A = 2,
	chr1B = 3,
	chr2 = 4,
	chr2A = 5,
	chr2a = 6,
	chr2B = 7,
	chr2b = 8,
	chr3 = 9,
	chr4 = 10,
	chr4A = 11,
	chr5 = 12,
	chr6 = 13,
	chr7 = 14,
	chr8 = 15,
	chr9 = 16,
	chr10 = 17,
	chr11 = 18,
	chr12 = 19,
	chr13 = 20,
	chr14 = 21,
	chr15 = 22,
	chr16 = 23,
	chr17 = 24,
	chr18 = 25,
	chr19 = 26,
	chr20 = 27,
	chr21 = 28,
	chr22 = 29,
	chr23 = 30,
	chr24 = 31,
	chr25 = 32,
	chr26 = 33,
	chr27 = 34,
	chr28 = 35,
	chr29 = 36,
	chr30 = 37,
	chr31 = 38,
	chr32 = 39,
	chr33 = 40,
	chr34 = 41,
	chr35 = 42,
	chr36 = 43,
	chr37 = 44,
	chr38 = 45,
	chr39 = 46,
	chr40 = 47,
	chrX = 48,
	chrY = 49,
	chrXY = 50,
	chrX1 = 51,
	chrX2 = 52,
	chrX3 = 53,
	chrX5 = 54,
	chrA1 = 55,
	chrA2 = 56,
	chrA3 = 57,
	chrB1 = 58,
	chrB2 = 59,
	chrB3 = 60,
	chrB4 = 61,
	chrC1 = 62,
	chrC2 = 63,
	chrD1 = 64,
	chrD2 = 65,
	chrD3 = 66,
	chrD4 = 67,
	chrE1 = 68,
	chrE2 = 69,
	chrE3 = 70,
	chrF1 = 71,
	chrF2 = 72,
	chrLG2 = 73,
	chrLG5 = 74,
	chrLGE22 = 75,
	chrW = 76,
	chrZ = 77,
	chrM = 78,
	chrLGE64 = 79,
	chrNA = -1,
	scaffold = 0
};

//possible strands.
enum Strand {
	fwd = 1,
	rev = -1,
	unk = 0
};

//possible frames
enum Frame {
	frame1 = 1,
	frame2 = 2,
	frame3 = 3,
	unknown = 0
};

//holds coordinates.
struct Coordinates {
	//the start position
	unsigned int start;
	//the end position
	unsigned int end;
	//nterm offset (see enum Offset)
	Offset Nterm;
	//cterm offset (see enum Offset)
	Offset Cterm;

	// lt opeator. returns true if lhs.start is lesser than rhs.start and lhs.end is smaller than rhs.end and lhs.end is smaller than rhs.start
	//otherwise returns false.
	bool operator()(const Coordinates& lhs, const Coordinates& rhs) const {
		return lhs.start < rhs.start && lhs.end < rhs.end && lhs.end < rhs.start;
	}
	//equal oparator. returns true if lhs.start is bigger than or equal to rhs.start and lhs.end smaller than or equal to rhs.end 
	//otherwise returns false.
	bool operator==(const Coordinates& rhs) const {
		return start >= rhs.start && end <= rhs.end;
	}
};

//extension of coordinates, holds genomic coordinates.
struct GenomeCoordinates : Coordinates {

	std::string transcriptid;
	std::string exonid;

	//holds the chromosome.
	Chromosome chr;
	//holds the scaffolding.
	std::string chrscaf;
	//holds the strand.
	Strand strand;
	//holds the frame.
	Frame frame;

	bool operator()(const GenomeCoordinates& lhs, const GenomeCoordinates& rhs) const {
		if (lhs.chr == scaffold && rhs.chr == scaffold && lhs.chrscaf == rhs.chrscaf) {
			return lhs.start < rhs.start && lhs.end < rhs.end && lhs.end >= rhs.start;
		}
		if (lhs.chr == scaffold && rhs.chr == scaffold && lhs.chrscaf != rhs.chrscaf) {
			return lhs.chrscaf < rhs.chrscaf;
		}
		if (lhs.chr == scaffold && rhs.chr != scaffold) {
			return false;
		}
		if (lhs.chr != scaffold && rhs.chr == scaffold) {
			return true;
		}
		if (lhs.chr == rhs.chr) {
			return lhs.start < rhs.start && lhs.end < rhs.end && lhs.end >= rhs.start;
		}
		return lhs.chr < rhs.chr;
	}

	bool operator==(const GenomeCoordinates& rhs) const {
		return ((chr != scaffold && chr == rhs.chr) || (chr == scaffold && chrscaf == rhs.chrscaf)) && start >= rhs.start && end <= rhs.end;
	}

	bool operator<(const GenomeCoordinates& rhs) const {
		if (chr == scaffold && rhs.chr == scaffold && chrscaf == rhs.chrscaf) {
			if (start == rhs.start) {
				return end < rhs.end;
			}
			return start < rhs.start;
		}
		if (chr == scaffold && rhs.chr == scaffold && chrscaf != rhs.chrscaf) {
			return chrscaf < rhs.chrscaf;
		}
		if (chr == scaffold && rhs.chr != scaffold) {
			return false;
		}
		if (chr != scaffold && rhs.chr == scaffold) {
			return true;
		}
		if (chr == rhs.chr) {
			if (start == rhs.start) {
				return end < rhs.end;
			}
			return start < rhs.start;
		}
		return chr < rhs.chr;
	}
};

#endif
