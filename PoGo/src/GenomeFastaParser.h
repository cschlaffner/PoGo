#ifndef GENOME_FASTA_PARSER_H
#define GENOME_FASTA_PARSER_H

#include "GTFParser.h"

class GenomeFastaParser {
public:
	static void readGenomeFASTA(std::string const& path);

};

#endif