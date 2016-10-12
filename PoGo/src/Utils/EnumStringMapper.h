#ifndef ENUMSTRINGMAPPER_H
#define ENUMSTRINGMAPPER_H

#include <iomanip>
#include <map>
#include "Coordinates.h"

//this static class contains methods to map the Enums in Utils.h 
//to strings or the other way around.
class EnumStringMapper {
public:
	static std::string enum_to_string(const Strand& strand, bool numeric = true);
	static std::string enum_to_string(const Chromosome& chr);
	static Strand string_to_strand(const std::string& string);
	static Frame string_to_frame(const std::string& string);
	static Chromosome string_to_chromosome(const std::string& string);
	static std::string ptm_to_colour(const std::string& ptmPSIname);
private:
	static std::map<std::string, std::string> m_ptm_to_colours;
};

#endif
