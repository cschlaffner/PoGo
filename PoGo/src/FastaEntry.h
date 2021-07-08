#ifndef FASTAENTRY_H
#define FASTAENTRY_H

#include "Utils.h"

class FastaEntry {
public:
	//ctr / dtr
	FastaEntry(std::string const& header, std::string const& AAsequence, std::string const& headersource);
	FastaEntry();
	~FastaEntry(void);
	//end ctr / dtr

	//returns the current header as string
	std::string get_header() const;
	//returns the current sequence.
	std::string get_sequence() const;
	//returns true if both header and sequence are empty
	bool is_empty() const;
	//returns the header source database
	std::string get_source() const;

private:
	//holds the fasta header
	std::string m_header;
	//holds the fasta entry sequence.
	std::string m_aa_sequence;
	//holds information about fasta header source
	std::string m_headersource;
};

#endif
