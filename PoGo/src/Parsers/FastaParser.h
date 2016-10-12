#ifndef FASTAPARSER_H
#define FASTAPARSER_H

#include <fstream>
//#include "FastaEntry.h"

//this fastaparser reads fasta files. 
//it will convert input sequences into iso-sequences [I, L] will be converted to J
class FastaParser {
public:
	//dtr
	~FastaParser();
	//end dtr.

	//singleton get_instance method.
	static FastaParser* get_instance();
	//opens the file
	bool open(const std::string& file);
	//closes the filestream
	void close();
	//parses and returns the next FASTA entry.
	FastaEntry next_entry();

	//overwritten assignment operator for singleton, so that no duplicates can be created.
	FastaParser& operator=(const FastaParser& other);
private:
	//inputstream.
	std::ifstream m_is;
	//current line.
	std::string m_line;
	//meyers singleton instance.
	static FastaParser* m_instance;
	//ctr
	FastaParser();
	FastaParser(FastaParser const&);
	//end ctr
};

class FastaParser__file_not_found_exception : public std::runtime_error {
public:
	FastaParser__file_not_found_exception(std::string error = "") : std::runtime_error(error) {};
};

#endif
