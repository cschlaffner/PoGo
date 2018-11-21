#include "GenomeFastaParser.h"

void GenomeFastaParser::readGenomeFASTA(std::string const& path) {
	std::ifstream	ifs(path.c_str());
	if (ifs.good()) {
		std::string line;
		while (std::getline(ifs, line)) {
			if (line.substr(0, 1) == ">") {
				std::size_t chrom_idx = line.find("dna:chromosome");
				std::size_t gscaf_idx = line.find("dna:chromosome");
				std::size_t scaff_idx = line.find("dna:chromosome");
				std::size_t white_idx = line.find(" ");
				if (chrom_idx != std::string::npos || gscaf_idx != std::string::npos) {
					if (white_idx != std::string::npos) {
						Chromosome::addChr(line.substr(1, white_idx - 1));
					}
				}
				else if (scaff_idx != std::string::npos) {
					if (white_idx != std::string::npos) {
						Chromosome::addScaffold(line.substr(1, white_idx - 1));
					}
				}
			}
		}
		GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::CHR_FROM_GENOME_FASTA = true;
		ifs.close();
	}
}