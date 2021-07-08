#include "FastaParser.h"

FastaParser* FastaParser::m_instance = nullptr;

FastaParser::FastaParser() { 
	m_headersource = "default";
	m_first_headercheck = true;
}

FastaParser::FastaParser(const FastaParser& other) {}
FastaParser::~FastaParser() {}

FastaParser* FastaParser::get_instance() {
	if (m_instance == nullptr) {
		m_instance = new FastaParser();
	}
	return m_instance;
}

bool FastaParser::open(const std::string& file) {
	if (!m_is.is_open()) {
		m_line = "";
		m_is.open(file.c_str());
	}
	return m_is.good();
}

void FastaParser::close() {
	m_line.clear();
	m_is.close();
	m_is.clear();
}

FastaEntry FastaParser::next_entry() {
	if (m_line.empty()) {
		std::getline(m_is, m_line);
	}
	std::string header = m_line;
	if (m_first_headercheck) {
		std::regex re(">(ENS[A-Z]*P.*)|(-)\\|(ENS[A-Z]*T.*)|(-)\\|(ENS[A-Z]*G.*)|(-)\\|(OTT[A-Z]*G.*)|(-)\\|(OTT[A-Z]*T.*)|(-)\\|.*");
		if (std::regex_match(header, re)==true) {
			m_headersource = "gencode";
			m_first_headercheck = false;
		} else {
			m_headersource = "other";
			m_first_headercheck = false;
		}
	}
	std::string sequence = "";
	while (std::getline(m_is, m_line) && (m_line.compare(0, 1, ">") != 0)) {
		sequence.append(make_iso_sequence(m_line));
	}
	if ((m_line.compare(0, 1, ">") != 0)) {
		m_line.clear();
	}
	return FastaEntry(header, sequence, m_headersource);
}

FastaParser& FastaParser::operator=(const FastaParser& other) {
	return *this;
}
