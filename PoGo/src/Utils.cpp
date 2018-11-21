#include "Utils.h"

std::string remove_ptms(const std::string& sequence) {
	std::string tmp = sequence;
	std::string::size_type start = tmp.find_first_of("(");
	while (start != std::string::npos) {
		std::string::size_type end = tmp.find_first_of(")");
		tmp.erase(tmp.begin() + start, tmp.begin() + end + 1);
		start = tmp.find_first_of("(");
	}

	std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

	return tmp;
}

std::string make_iso_sequence(std::string sequence) {
	std::replace(sequence.begin(), sequence.end(), 'I', 'J');
	std::replace(sequence.begin(), sequence.end(), 'L', 'J');

	return sequence;
}

void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters, bool trimEmpty) {
	typedef std::vector<std::string> container_t;
	typedef container_t::value_type value_type;
	typedef container_t::size_type size_type;

	std::string::size_type pos, last_pos = 0;

	while (true) {
		pos = str.find_first_of(delimiters, last_pos);
		if (pos == std::string::npos) {
			pos = str.length();

			if (pos != last_pos || !trimEmpty) {
				tokens.push_back(value_type(str.data() + last_pos, pos - last_pos));
			}
			break;
		}
		if (pos != last_pos || !trimEmpty) {
			tokens.push_back(value_type(str.data() + last_pos, size_type(pos) - last_pos));
		}
		last_pos = pos + 1;
	}
}

CoordinateMapType::value_type get_coordinates(Coordinates const& proteinCoords, GenomeCoordinates& genomeCoords, Coordinates& peptideCoords) {
	unsigned int pep_start;
	unsigned int pep_end;
	
	if (peptideCoords.start >= proteinCoords.start && peptideCoords.start <= proteinCoords.end) {
		//if the start of the peptide lies within the protein use this as the protein start
		pep_start = peptideCoords.start;
	} else {
		//otherwise the peptide has to start where the protein starts.
		pep_start = proteinCoords.start;
	}

	if (peptideCoords.end >= proteinCoords.start && peptideCoords.end <= proteinCoords.end) {
		//if the end of the also lies within the protein use this as peptide end
		pep_end = peptideCoords.end;
	} else {
		//otherwise the peptide has to end where the protein ends
		pep_end = proteinCoords.end;
	}
	Coordinates partial_peptide_coords;
	partial_peptide_coords.start = pep_start;
	partial_peptide_coords.end = pep_end;
	int start_genomic_coord;
	int end_genomic_coord;

	if (genomeCoords.end >= genomeCoords.start) {
		//coding from forward strand
		if (pep_start == proteinCoords.start) {
			//peptide starting at protein start position
			partial_peptide_coords.Nterm = proteinCoords.Nterm;
			start_genomic_coord = genomeCoords.start;
		} else {
			//peptide starting within the protein
			partial_peptide_coords.Nterm = Offset(3);
			start_genomic_coord = ((pep_start - 1 - proteinCoords.start) * 3) + genomeCoords.start + int(proteinCoords.Nterm);
		}
		if (pep_end == proteinCoords.end) {
			//peptide ending at protein end position
			partial_peptide_coords.Cterm = proteinCoords.Cterm;
			end_genomic_coord = genomeCoords.end;
		} else {
			//peptide ending within protein
			partial_peptide_coords.Cterm = Offset(3);
			end_genomic_coord = ((pep_end - proteinCoords.start) * 3) + genomeCoords.start + (int(proteinCoords.Nterm) - 1);
		}
	}
	else {
		// coding from reverse strand! 
		if (pep_start == proteinCoords.start) {
			//peptide starting at protein start position
			partial_peptide_coords.Nterm = proteinCoords.Nterm;
			start_genomic_coord = genomeCoords.start;
		} else {
			//peptide starting within the protein
			partial_peptide_coords.Nterm = Offset(3);
			start_genomic_coord = genomeCoords.start - (((pep_start - proteinCoords.start - 1) * 3) + int(proteinCoords.Nterm));
		}
		if (pep_end == proteinCoords.end) {
			//peptide ending at protein end position
			partial_peptide_coords.Cterm = proteinCoords.Cterm;
			end_genomic_coord = genomeCoords.end;
		} else {
			//peptide ending within protein
			partial_peptide_coords.Cterm = Offset(3);
			end_genomic_coord = genomeCoords.start - (((pep_end - proteinCoords.start) * 3) + int(proteinCoords.Nterm) - 1);
		}
	}

	GenomeCoordinates genome_coordinates;
	genome_coordinates.chr = genomeCoords.chr;
	genome_coordinates.chrscaf = genomeCoords.chrscaf;
	genome_coordinates.strand = genomeCoords.strand;
	genome_coordinates.frame = genomeCoords.frame;
	genome_coordinates.start = start_genomic_coord;
	genome_coordinates.end = end_genomic_coord;
	genome_coordinates.transcriptid = genomeCoords.transcriptid;
	genome_coordinates.exonid = genomeCoords.exonid;

	//flip start and end if coding from reverse strand
	if (genome_coordinates.strand == rev) {
		unsigned int start = genome_coordinates.end;
		unsigned int end = genome_coordinates.start;
		genome_coordinates.start = start;
		genome_coordinates.end = end;
	}

	return CoordinateMapType::value_type(partial_peptide_coords, genome_coordinates);
}

bool same_coordinates(const GenomeCoordinates& lhs, const GenomeCoordinates& rhs) {
	return ((!lhs.chr.isScaffold() && lhs.chr.getValue() == rhs.chr.getValue())
		|| (lhs.chr.isScaffold() && lhs.chrscaf == rhs.chrscaf))
		&& (lhs.start == rhs.start)
		&& (lhs.end == rhs.end)
		&& (lhs.frame == rhs.frame)
		&& (lhs.strand == rhs.strand);
}

bool compare_coordinates_ascending(const GenomeCoordinates& lhs, const GenomeCoordinates& rhs) {
	return lhs.start < rhs.start;
}

std::string coordinates_to_string(const GenomeCoordinates& coords, bool chrincluded) {
	std::stringstream ss;
	if (coords.chr.isScaffold() && !chrincluded) {
		ss << coords.chrscaf;
	}else if(coords.chr.isScaffold() && chrincluded){
		ss << "scaffold" << coords.chrscaf;
	} else if(!chrincluded){
		ss << EnumStringMapper::enum_to_string(coords.chr);
	} else if(chrincluded){
		ss << EnumStringMapper::enum_to_chr_string(coords.chr);
	}
	ss << ":" << coords.start << "-" << coords.end << " " << EnumStringMapper::enum_to_string(coords.strand, false);// << " " << coords.frame;
	return ss.str();
}

std::string coordinates_to_short_string(const GenomeCoordinates& coords, unsigned int offset, bool chrincluded) {
	std::stringstream ss;
	if (coords.chr.isScaffold() && !chrincluded) {
		ss << coords.chrscaf;
	}else if(coords.chr.isScaffold() && chrincluded){
		ss << "scaffold" << coords.chrscaf;
	} else if(!chrincluded){
		ss << EnumStringMapper::enum_to_string(coords.chr);
	} else if(chrincluded){
		ss << EnumStringMapper::enum_to_chr_string(coords.chr);
	}
	ss << ":" << (coords.start - offset) << "-" << coords.end;
	return ss.str();
}

std::string coordinates_to_gtf_string(const GenomeCoordinates& coords, const std::string& type, bool frameinclude, const std::string& source, bool chrincluded) {
	std::stringstream ss;
	if (coords.chr.isScaffold() && !chrincluded) {
		ss << coords.chrscaf;
	}else if(coords.chr.isScaffold() && chrincluded){
		ss << "scaffold" << coords.chrscaf;
	} else if(!chrincluded){
		ss << EnumStringMapper::enum_to_string(coords.chr);
	} else if(chrincluded){
		ss << EnumStringMapper::enum_to_chr_string(coords.chr);
	}
	ss << "\t" << source << "\t" << type << "\t" << coords.start << "\t" << coords.end << "\t.\t" << EnumStringMapper::enum_to_string(coords.strand, false);
	if (frameinclude) {
		ss << "\t" << coords.frame << "\t";
	} else {
		ss << "\t.\t";
	}
	return ss.str();
}

std::string coordinates_to_bed_string(const GenomeCoordinates& coords, const std::string& name, unsigned int score, bool chrincluded) {
	std::stringstream ss;
	if (coords.chr.isScaffold() && !chrincluded) {
		ss << coords.chrscaf;
	}else if(coords.chr.isScaffold() && chrincluded){
		ss << "scaffold" << coords.chrscaf;
	} else if(!chrincluded){
		ss << EnumStringMapper::enum_to_string(coords.chr);
	} else if(chrincluded){
		ss << EnumStringMapper::enum_to_chr_string(coords.chr);
	}
	ss << "\t" << (coords.start - 1) << "\t" << coords.end << "\t" << name << "\t" << score << "\t" << EnumStringMapper::enum_to_string(coords.strand, false) << "\t" << (coords.start - 1) << "\t" << (coords.start - 1) << "\t";
	return ss.str();
}

std::string coordinates_to_short_bed_string(const GenomeCoordinates& coords, const std::string& name, unsigned int score, bool chrincluded) {
	std::stringstream ss;
	if (coords.chr.isScaffold() && !chrincluded) {
		ss << coords.chrscaf;
	}else if(coords.chr.isScaffold() && chrincluded){
		ss << "scaffold" << coords.chrscaf;
	} else if(!chrincluded){
		ss << EnumStringMapper::enum_to_string(coords.chr);
	} else if(chrincluded){
		ss << EnumStringMapper::enum_to_chr_string(coords.chr);
	}

	ss << "\t" << (coords.start - 1) << "\t" << coords.end << "\t" << name << "\t" << score << "\t" << EnumStringMapper::enum_to_string(coords.strand, false) << "\t";
	return ss.str();
}

std::string coordinates_to_gct_string(std::vector<GenomeCoordinates> const& coords, bool chrincluded) {
	std::stringstream ss;
	for (size_t i(0); i < coords.size(); ++i) {
		if (i > 0) {
			ss << ",";
		}
		ss << coordinates_to_short_string(coords.at(i), 1, chrincluded);
	}
	return ss.str();
}

GenomeCoordinates extract_coordinates_from_gtf_line(std::vector<std::string> const&tokens) {
	GenomeCoordinates coord;
	coord.chr = EnumStringMapper::string_to_chromosome(tokens.at(0));
	if (coord.chr.isScaffold()) {
		coord.chrscaf = tokens.at(0);
	} else {
		coord.chrscaf = "";
	}
	coord.strand = EnumStringMapper::string_to_strand(tokens.at(6));
	coord.frame = EnumStringMapper::string_to_frame(tokens.at(7));
	if (coord.strand == fwd) {
		coord.start = atoi(tokens.at(3).c_str());
		coord.end = atoi(tokens.at(4).c_str());
	} else if (coord.strand == rev) {
		coord.end = atoi(tokens.at(3).c_str());
		coord.start = atoi(tokens.at(4).c_str());
	}
	return coord;
}

bool compare_genome_coordinate_sets_ascending(const std::vector<GenomeCoordinates>& lhs, const std::vector<GenomeCoordinates>& rhs) {
	bool same(false);
	size_t lhs_it(0);
	size_t rhs_it(0);
	if (lhs.size() > 1 && (lhs.size() == rhs.size())) {
		//if both sets have the same amount of entries
		while (lhs_it < lhs.size() && rhs_it < rhs.size() && (!same)) {
			//compare the entries. 
			same = compare_coordinates_ascending_whole(lhs.at(lhs_it), rhs.at(rhs_it));
			++lhs_it;
			++rhs_it;
		}
	} else if (lhs.size() == 1 && rhs.size() == 1) {
		//if both sets only contain one entry
		same = compare_coordinates_ascending_whole(lhs.at(lhs_it), rhs.at(rhs_it));
	} else {
		//if one contains less entries than another
		size_t it_end = (lhs.size() < rhs.size()) ? lhs.size() : rhs.size();
		for (size_t i(0); (i < it_end) && (!same); ++i) {
			same = compare_coordinates_ascending_whole(lhs.at(i), rhs.at(i));
		}
	}
	return same;
}

bool compare_coordinates_ascending_whole(const GenomeCoordinates& lhs, const GenomeCoordinates& rhs) {
	return lhs.start < rhs.start || (lhs.start == rhs.start && lhs.end < rhs.end);
}

bool byIntValue::operator()(const std::pair<std::string, unsigned>& lhs, const std::pair<std::string, unsigned>& rhs) const {
	return lhs.second <= rhs.second;
}

/**
 * This function check that the filename ends on the String name.
 * @param nameFile name of the file
 * @param extension extension of the file
 * @return return true if the exntension is the externsion of the file.
 */
bool isInLastPosition(std::string nameFile, std::string extension){
	return nameFile.rfind(extension) == (nameFile.size() - extension.size());
}

/**
 * This function remove the extension of the file if is found, if the extension is not found, it will keep
 * the original value.
 * @param nameFile Name of the file.
 * @param extension Extension to be removed
 * @return new File Name
 */
std::string removeExtensionOutput(std::string nameFile, std::string extension){
	std::size_t found = nameFile.rfind(extension);
	if (found!=std::string::npos)
		nameFile.replace(found,extension.length(),"");

	return nameFile;
}