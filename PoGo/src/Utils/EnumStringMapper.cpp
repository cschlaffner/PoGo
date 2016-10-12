//#include "EnumStringMapper.h"

std::map<std::string, std::string> EnumStringMapper::m_ptm_to_colours(
	std::map<std::string, std::string>({ 
		std::pair<std::string,  std::string>("phospho",		"255,51,51"),
		std::pair<std::string,  std::string>("acetyl",		"204,102,0"),
		std::pair<std::string,  std::string>("amidated",	"255,153,51"),
		std::pair<std::string,  std::string>("oxidation",	"204,204,0"),
		std::pair<std::string,  std::string>("methyl",		"0,204,0"),
		std::pair<std::string,  std::string>("glygly",		"51,255,51"),
		std::pair<std::string,  std::string>("gg",			"51,255,51"),
		std::pair<std::string,  std::string>("sulfo",		"0,204,204"),
		std::pair<std::string,  std::string>("palmitoyl",	"51,153,255"),
		std::pair<std::string,  std::string>("formyl",		"0,0,204"),
		std::pair<std::string,  std::string>("deamidated",  "51,51,255") 
	}));

std::string EnumStringMapper::enum_to_string(const Strand& strand, bool numeric) {
	if (numeric == true) {
		switch (strand) {
		case fwd: return "1";
		case rev: return "-1";
		default: return "0";
		}
	}
	switch (strand) {
	case fwd: return "+";
	case rev: return "-";
	default: return ".";
	}
}

std::string EnumStringMapper::enum_to_string(const Chromosome& chr) {
	switch (chr) {
	case chr1: return "chr1";
	case chr2: return "chr2";
	case chr3: return "chr3";
	case chr4: return "chr4";
	case chr5: return "chr5";
	case chr6: return "chr6";
	case chr7: return "chr7";
	case chr8: return "chr8";
	case chr9: return "chr9";
	case chr10: return "chr10";
	case chr11: return "chr11";
	case chr12: return "chr12";
	case chr13: return "chr13";
	case chr14: return "chr14";
	case chr15: return "chr15";
	case chr16: return "chr16";
	case chr17: return "chr17";
	case chr18: return "chr18";
	case chr19: return "chr19";
	case chr20: return "chr20";
	case chr21: return "chr21";
	case chr22: return "chr22";
	case chrX: return "chrX";
	case chrY: return "chrY";
	case chrXY: return "chrXY";
	case chrM: return "chrM";
	case scaffold: return "";
	default: return "unknown";
	}
}

Strand EnumStringMapper::string_to_strand(const std::string& string) {
	if (string.compare("-1") == 0 || string.compare("-") == 0) {
		return rev;
	}
	if (string.compare("1") == 0 || string.compare("+") == 0) {
		return fwd;
	}
	return unk;
}

Frame EnumStringMapper::string_to_frame(const std::string& string) {
	if (string.compare("1") == 0) {
		return frame1;
	}
	if (string.compare("2") == 0) {
		return frame2;
	}
	if (string.compare("3") == 0 || string.compare("0") == 0) {
		return frame3;
	}
	return unknown;
}

Chromosome EnumStringMapper::string_to_chromosome(const std::string& string) {
	std::string substr;
	if (string.size() > 3 && (string.compare(0, 3, "chr") == 0 || string.compare(0, 3, "Chr") == 0)) {
		substr = string.substr(3, string.size() - 1);
	} else {
		substr = string;
	}

	if (substr.compare("X") == 0) {
		return chrX;
	}
	if (substr.compare("Y") == 0) {
		return chrY;
	}
	if (substr.compare("XY") == 0) {
		return chrXY;
	}
	if (substr.compare("M") == 0 || substr.compare("MT") == 0) {
		return chrM;
	}
	if (substr.substr(0, 2).compare("GL") == 0 || substr.substr(0, 2).compare("KI") == 0 || substr.substr(0, 2).compare("JH") == 0 || substr.substr(0, 2).compare("KB") == 0) {
		return scaffold;
	}
	int chrom_num = atoi(substr.c_str());
	return Chromosome(chrom_num);
}

std::string EnumStringMapper::ptm_to_colour(const std::string& ptmPSIname) {
	std::string name = ptmPSIname;
	for (size_t i(0); i < name.size(); i++) {
		name[i] = tolower(name[i]);
	}
	if (m_ptm_to_colours.count(name) == 0) {
		return "255,51,153";
	}
	return m_ptm_to_colours.at(name);
}
