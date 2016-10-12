#ifndef COORDINATES_H
#define COORDINATES_H

//#include "Globals.h"

//possible offsets.
enum Offset {
	off1 = 1,
	off2 = 2,
	off3 = 3
};

//possible chromosomes
enum Chromosome {
	chr1 = 1,
	chr2 = 2,
	chr3 = 3,
	chr4 = 4,
	chr5 = 5,
	chr6 = 6,
	chr7 = 7,
	chr8 = 8,
	chr9 = 9,
	chr10 = 10,
	chr11 = 11,
	chr12 = 12,
	chr13 = 13,
	chr14 = 14,
	chr15 = 15,
	chr16 = 16,
	chr17 = 17,
	chr18 = 18,
	chr19 = 19,
	chr20 = 20,
	chr21 = 21,
	chr22 = 22,
	chrX = 23,
	chrY = 24,
	chrXY = 25,
	chrM = 26,
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
