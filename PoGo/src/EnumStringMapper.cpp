#include "EnumStringMapper.h"

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
	if (numeric) {
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
	std::string name = chr.getName();
	if (name.size() > 3) {
		if (name.substr(0, 3) == "chr" || name.substr(0, 3) == "Chr") {
			return name.substr(3);
		}
	}
	return name;
}

std::string EnumStringMapper::enum_to_chr_string(const Chromosome& chr) {
	std::string name = enum_to_string(chr);
	if (name.size() < 3) {
		return "chr" + name;
	}
	else {
		if (name.substr(0, 3) != "chr" && name.substr(0, 3) != "Chr" && !chr.isScaffold() && !chr.isNA()) {
			return "chr" + name;
		}
	}
	return name;
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
	}
	else {
		substr = string;
	}

	if (GENOME_MAPPER_GLOBALS::PEPTIDE_MAPPER::CHR_FROM_GENOME_FASTA == false) {
		Chromosome::addChr(substr);
	}
	return Chromosome(substr);
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
