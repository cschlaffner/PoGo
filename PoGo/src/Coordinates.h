#ifndef COORDINATES_H
#define COORDINATES_H

#include "Chromosome.h"

//possible offsets.
enum Offset {
	off1 = 1,
	off2 = 2,
	off3 = 3
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

	Coordinates() {
		start = 0;
		end = 0;
		Nterm = Offset::off3;
		Cterm = Offset::off3;
	}

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

	GenomeCoordinates() : Coordinates() {
		chr = Chromosome();
		chrscaf = "";
		strand = Strand::unk;
		frame = Frame::unknown;
	}

	bool operator()(const GenomeCoordinates& lhs, const GenomeCoordinates& rhs) const {
		if (lhs.chr.isScaffold() && rhs.chr.isScaffold() && lhs.chrscaf == rhs.chrscaf) {
			return lhs.start < rhs.start && lhs.end < rhs.end && lhs.end >= rhs.start;
		}
		if (lhs.chr.isScaffold() && rhs.chr.isScaffold() && lhs.chrscaf != rhs.chrscaf) {
			return lhs.chrscaf < rhs.chrscaf;
		}
		if (lhs.chr.isScaffold() && !rhs.chr.isScaffold()) {
			return false;
		}
		if (!lhs.chr.isScaffold() && rhs.chr.isScaffold()) {
			return true;
		}
		if (lhs.chr.getValue() == rhs.chr.getValue()) {
			if ((lhs.strand == Strand::rev && rhs.strand == Strand::rev) || (lhs.strand == Strand::unk && rhs.strand == Strand::rev) || (lhs.strand == Strand::rev && rhs.strand == Strand::unk)) {
				return lhs.start > rhs.start && lhs.end > rhs.end && lhs.start >= rhs.end;
			} else {
				return lhs.start < rhs.start && lhs.end < rhs.end && lhs.end <= rhs.start;
			}
		}
		return lhs.chr.getValue() < rhs.chr.getValue();
	}

	bool operator==(const GenomeCoordinates& rhs) const {
		return ((!chr.isScaffold() && chr.getValue() == rhs.chr.getValue()) || (chr.isScaffold() && chrscaf == rhs.chrscaf)) && start >= rhs.start && end <= rhs.end;
	}

	bool operator<(const GenomeCoordinates& rhs) const {
		if (chr.isScaffold() && rhs.chr.isScaffold() && chrscaf == rhs.chrscaf) {
			if ((strand == Strand::rev && rhs.strand == Strand::rev) || (strand == Strand::unk && rhs.strand == Strand::rev) || (strand == Strand::rev && rhs.strand == Strand::unk)) {
				if (end == rhs.end) {
					return start > rhs.start;
				}
				return end > rhs.end;
			} else {
				if (start == rhs.start) {
					return end < rhs.end;
				}
				return start < rhs.start;
			}
		}
		if (chr.isScaffold() && rhs.chr.isScaffold() && chrscaf != rhs.chrscaf) {
			return chrscaf < rhs.chrscaf;
		}
		if (chr.isScaffold() && !rhs.chr.isScaffold()) {
			return false;
		}
		if (!chr.isScaffold() && rhs.chr.isScaffold()) {
			return true;
		}
		if (chr.getValue() == rhs.chr.getValue()) {
			if ((strand == Strand::rev && rhs.strand == Strand::rev) || (strand == Strand::unk && rhs.strand == Strand::rev) || (strand == Strand::rev && rhs.strand == Strand::unk)) {
				if (end == rhs.end) {
					return start > rhs.start;
				}
				return end > rhs.end;
			} else {
				if (start == rhs.start) {
					return end < rhs.end;
				}
				return start < rhs.start;
			}
		}
		return chr.getValue() < rhs.chr.getValue();
	}
};

#endif
