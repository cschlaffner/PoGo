#include "GeneEntry.h"

//PUBLIC
//ctr & dtr
GeneEntry::GeneEntry(void) {}
GeneEntry::GeneEntry(std::string const& gtfgeneline) {
	init(gtfgeneline);
}
GeneEntry::~GeneEntry(void) {}

std::string GeneEntry::get_id() const {
	return m_id;
}

std::string GeneEntry::get_type() const {
	return m_type;
}

std::string GeneEntry::get_status() const {
	return m_status;
}

std::string GeneEntry::get_name() const {
	return m_gene_name;
}

bool GeneEntry::operator<(const GeneEntry& rhs) const {
	if (to_string(int(m_coord.chr)) < to_string(int(rhs.m_coord.chr))) {
		return true;
	}
	if (m_coord.chr == rhs.m_coord.chr) {
		if (m_coord.start == rhs.m_coord.start) {
			return m_coord.end < rhs.m_coord.end;
		}
		return m_coord.start < rhs.m_coord.start;
	}
	return false;
}

std::ostream& GeneEntry::to_gtf(const std::string& source, std::ostream& os, bool chrincluded) {
	os << coordinates_to_gtf_string(m_coord, "gene", false, source, chrincluded);
	os << "gene_id \"" << m_id << "\"; transcript_id \"" << m_id << "\"; gene_type \"" << m_type << "\"; gene_status \"" << m_status << "\"; gene_name \"" << m_gene_name;
	os << "\"; transcript_type \"" << m_type << "\"; transcript_status \"" << m_status << "\"; transcript_name \"" << m_gene_name << "\";";

	for (size_t i(0); i < m_tags.size(); ++i) {
		os << " tag \"" << m_tags.at(i) << "\";";
	}
	return os;
}

bool GeneEntry::is_primary() {
	return m_coord.chr != chrNA && m_coord.chr != scaffold&& m_coord.chr != chrXY;
}

bool GeneEntry::is_patchhaploscaff() {
	return m_coord.chr == scaffold && m_coord.chrscaf != "";
}

std::string GeneEntry::extract_gene_id(std::string gtfGeneLine) {
	std::size_t index = gtfGeneLine.find(GENOME_MAPPER_GLOBALS::ID::GTF_GENE_ID) + GENOME_MAPPER_GLOBALS::ID::GTF_GENE_ID.length();
	std::string value = "";
	if (index != std::string::npos) {
		while (gtfGeneLine[index] != '\"') {
			value = value + gtfGeneLine[index];
			index += 1;
		}
	}
	return value;
}

std::string GeneEntry::extract_transcript_id(std::string gtfGeneLine) {
	std::size_t index = gtfGeneLine.find(GENOME_MAPPER_GLOBALS::ID::GTF_TRANSCRIPT_ID) + GENOME_MAPPER_GLOBALS::ID::GTF_TRANSCRIPT_ID.length();
	std::string value = "";
	if (index != std::string::npos) {
		while (gtfGeneLine[index] != '\"') {
			value = value + gtfGeneLine[index];
			index += 1;
		}
	}
	return value;
}

std::string GeneEntry::extract_exon_id(std::string gtfGeneLine) {
	std::size_t index = gtfGeneLine.find(GENOME_MAPPER_GLOBALS::ID::GTF_EXON_ID) + GENOME_MAPPER_GLOBALS::ID::GTF_EXON_ID.length();
	std::string value = "";
	if (index != std::string::npos) {
		if (gtfGeneLine[index] != '\"') {
			value = value + gtfGeneLine[index];
			index += 1;
		}
	}
	return value;
}

//PRIVATE
void GeneEntry::init(const std::string& gtfgeneline) {
	std::vector<std::string> tokens;
	tokenize(gtfgeneline, tokens, "\t");

	init(extract_gene_id(gtfgeneline), extract_coordinates_from_gtf_line(tokens), extract_type(tokens), extract_status(tokens), extract_gene_name(tokens), extract_tags(tokens));
}

void GeneEntry::init(std::string ID, GenomeCoordinates coordinates, std::string type, std::string status, std::string gene_name, std::vector<std::string> tags) {
	m_coord = coordinates;
	m_id = ID;
	m_type = type;
	m_status = status;
	m_gene_name = gene_name;
	m_tags = tags;
}

std::string GeneEntry::extract_type(const std::vector<std::string>& tokens) {
	std::string value;
	if (tokens.size() >= 9) {
		std::vector<std::string> res = extract_by_tag("gene_type", tokens.at(8));
		if (res.size() == 1) {
			value = res.at(0);
		}
	}
	return value;
}

std::string GeneEntry::extract_status(const std::vector<std::string>& tokens) {
	std::string value;
	if (tokens.size() >= 9) {
		std::vector<std::string> res = extract_by_tag("gene_status", tokens.at(8));
		if (res.size() == 1) {
			value = res.at(0);
		}
	}
	return value;
}

std::string GeneEntry::extract_gene_name(const std::vector<std::string>& tokens) {
	std::string value;
	if (tokens.size() >= 9) {
		std::vector<std::string> res = extract_by_tag("gene_name", tokens.at(8));
		if (res.size() == 1) {
			value = res.at(0);
		}
	}
	return value;
}

std::vector<std::string> GeneEntry::extract_tags(const std::vector<std::string>& tokens) {
	std::vector<std::string> values;
	if (tokens.size() >= 9) {
		values = extract_by_tag("tag", tokens.at(8));
	}
	return values;
}

std::vector<std::string> GeneEntry::extract_by_tag(const std::string& tag, const std::string& tagList) {
	std::vector<std::string> values = std::vector<std::string>();
	std::vector<std::string> return_values = std::vector<std::string>();

	tokenize(tagList, values, ";", true);
	for (size_t i(0); i < values.size(); ++i) {
		int start = 0;
		if (values.at(i).substr(0, 1).compare(" ") == 0) {
			start = 1;
		}
		if (values.at(i).substr(start, tag.size()).compare(tag) == 0) {
			std::vector<std::string> v = std::vector<std::string>();
			tokenize(values.at(i), v, "\"");
			if (v.size() >= 2) {
				return_values.push_back(v.at(1));
			}
		}
	}
	return return_values;
}

