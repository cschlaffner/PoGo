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
	switch (chr) {
	case chr1: return "1";
	case chr1A: return "1A";
	case chr1B: return "1B";
	case chr2: return "2";
	case chr2A: return "2A";
	case chr2a: return "2a";
	case chr2B: return "2B";
	case chr2b: return "2b";
	case chr3: return "3";
	case chr4: return "4";
	case chr4A: return "4A";
	case chr5: return "5";
	case chr6: return "6";
	case chr7: return "7";
	case chr8: return "8";
	case chr9: return "9";
	case chr10: return "10";
	case chr11: return "11";
	case chr12: return "12";
	case chr13: return "13";
	case chr14: return "14";
	case chr15: return "15";
	case chr16: return "16";
	case chr17: return "17";
	case chr18: return "18";
	case chr19: return "19";
	case chr20: return "20";
	case chr21: return "21";
	case chr22: return "22";
	case chr23: return "23";
	case chr24: return "24";
	case chr25: return "25";
	case chr26: return "26";
	case chr27: return "27";
	case chr28: return "28";
	case chr29: return "29";
	case chr30: return "30";
	case chr31: return "31";
	case chr32: return "32";
	case chr33: return "33";
	case chr34: return "34";
	case chr35: return "35";
	case chr36: return "36";
	case chr37: return "37";
	case chr38: return "38";
	case chr39: return "39";
	case chr40: return "40";
	case chrX: return "chrX";
	case chrY: return "chrY";
	case chrXY: return "chrXY";
	case chrX1: return "X1";
	case chrX2: return "X2";
	case chrX3: return "X3";
	case chrX5: return "X5";
	case chrA1: return "A1";
	case chrA2: return "A2";
	case chrA3: return "A3";
	case chrB1: return "B1";
	case chrB2: return "B2";
	case chrB3: return "B3";
	case chrB4: return "B4";
	case chrC1: return "C1";
	case chrC2: return "C2";
	case chrD1: return "D1";
	case chrD2: return "D2";
	case chrD3: return "D3";
	case chrD4: return "D4";
	case chrE1: return "E1";
	case chrE2: return "E2";
	case chrE3: return "E3";
	case chrF1: return "F1";
	case chrF2: return "F2";
	case chrLG2: return "LG2";
	case chrLG5: return "LG5";
	case chrLGE22: return "LGE22";
	case chrW: return "W";
	case chrZ: return "Z";
	case chrM: return "MT";
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
	if (substr.compare("W") == 0) {
		return chrW;
	}
	if (substr.compare("Z") == 0) {
		return chrZ;
	}
	if (substr.compare("A1") == 0) {
		return chrA1;
	}
	if (substr.compare("A2") == 0) {
		return chrA2;
	}
	if (substr.compare("A3") == 0) {
		return chrA3;
	}
	if (substr.compare("B1") == 0) {
		return chrB1;
	}
	if (substr.compare("B2") == 0) {
		return chrB2;
	}
	if (substr.compare("B3") == 0) {
		return chrB3;
	}
	if (substr.compare("B4") == 0) {
		return chrB4;
	}
	if (substr.compare("C1") == 0) {
		return chrC1;
	}
	if (substr.compare("C2") == 0) {
		return chrC2;
	}
	if (substr.compare("D1") == 0) {
		return chrD1;
	}
	if (substr.compare("D2") == 0) {
		return chrD2;
	}
	if (substr.compare("D3") == 0) {
		return chrD3;
	}
	if (substr.compare("D4") == 0) {
		return chrD4;
	}
	if (substr.compare("E1") == 0) {
		return chrE1;
	}
	if (substr.compare("E2") == 0) {
		return chrE2;
	}
	if (substr.compare("E3") == 0) {
		return chrE3;
	}
	if (substr.compare("F1") == 0) {
		return chrF1;
	}
	if (substr.compare("F2") == 0) {
		return chrF2;
	}
	if (substr.compare("LGE64") == 0) {
		return chrLGE64;
	}
	if (substr.compare("2A") == 0) {
		return chr2A;
	}
	if (substr.compare("2B") == 0) {
		return chr2B;
	}
	if (substr.compare("2a") == 0) {
		return chr2a;
	}
	if (substr.compare("2b") == 0) {
		return chr2b;
	}
	if (substr.compare("X1") == 0) {
		return chrX1;
	}
	if (substr.compare("X2") == 0) {
		return chrX2;
	}
	if (substr.compare("X3") == 0) {
		return chrX3;
	}
	if (substr.compare("X5") == 0) {
		return chrX5;
	}
	if (substr.compare("1A") == 0) {
		return chr1A;
	}
	if (substr.compare("1B") == 0) {
		return chr1B;
	}
	if (substr.compare("4A") == 0) {
		return chr4A;
	}
	if (substr.compare("LG2") == 0) {
		return chrLG2;
	}
	if (substr.compare("LG5") == 0) {
		return chrLG5;
	}
	if (substr.compare("LGE22") == 0) {
		return chrLGE22;
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
	if (substr.length() >= 2) {
		if (substr.substr(0, 2).compare("GL") == 0 || substr.substr(0, 2).compare("KI") == 0 || substr.substr(0, 2).compare("JH") == 0 || substr.substr(0, 2).compare("KB") == 0 || substr.substr(0, 2).compare("GJ") == 0 || substr.substr(0, 2).compare("HT") == 0 || substr.substr(0, 2).compare("KL") == 0 || substr.substr(0, 2).compare("KE") == 0 || substr.substr(0, 2).compare("Un") == 0 || substr.substr(0, 2).compare("un") == 0 || substr.substr(0, 2).compare("cu") == 0 || substr.substr(0, 2).compare("KQ") == 0 || substr.substr(0, 2).compare("sc") == 0 || substr.substr(0, 2).compare("ul") == 0 || substr.substr(0, 2).compare("Co") == 0 || substr.substr(0, 2).compare("Ul") == 0) {
			return scaffold;
		}
	}
	if (substr.length() >= 4) {
		if (substr.substr(0, 4).compare("ACFV") == 0 || substr.substr(0, 4).compare("AAEX") == 0 || substr.substr(0, 4).compare("AQIB") == 0 || substr.substr(0, 4).compare("JSUE") == 0 || substr.substr(0, 4).compare("AAGW") == 0 || substr.substr(0, 4).compare("AMGL") == 0 || substr.substr(0, 4).compare("AACZ") == 0 || substr.substr(0, 4).compare("AHZZ") == 0 || substr.substr(0, 4).compare("AABR") == 0 || substr.substr(0, 4).compare("AAND") == 0 || substr.substr(substr.size() - 4, 4).compare("hap1") == 0 || substr.substr(substr.size() - 4, 4).compare("hap2") == 0 || substr.substr(substr.size() - 4, 4).compare("ndom") == 0) {
			return scaffold;
		}
	}
	int chrom_num = atoi(substr.c_str());
	if (chrom_num == 2) {
		chrom_num = 4;
	}
	else if (chrom_num == 3 || chrom_num == 4) {
		chrom_num = chrom_num + 6;
	}
	else {
		chrom_num = chrom_num + 7;
	}
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
