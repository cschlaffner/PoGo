#ifndef RESULTPARSER_H
#define RESULTPARSER_H

//#include "Utils\MappedPeptides.h"

//this class parses peptide input files and maps them to the genome.
class ResultParser {
public:
	//dtr
	~ResultParser(void);
	//end dtr

	//read function. this reads the peptides input and sets the wheels in motion.
	void static read(std::string file, CoordinateWrapper& coordwrapper, MappedPeptides& mapping, std::string unmappedoutput, KmereMap& k);
private:
	//ctr
	ResultParser(void);
	// end ctr.
};

class ResultParser__output_file_exception : public std::runtime_error {
public:
	ResultParser__output_file_exception(std::string err = "") : std::runtime_error(err) {}
};



#endif
