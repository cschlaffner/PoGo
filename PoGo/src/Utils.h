#ifndef UTILS_H
#define UTILS_H

#include <set>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include "EnumStringMapper.h"
//---------------------------
//Utils.h and Utils.cpp hold an abundance of useful functions.
//---------------------------
typedef std::multimap<Coordinates, GenomeCoordinates, Coordinates> CoordinateMapType;
//this function creates the an entry for a CoordinateMapType.
CoordinateMapType::value_type get_coordinates(Coordinates const& proteinCoords, GenomeCoordinates& genomeCoords, Coordinates& peptideCoords);

//removes all ptms from a sequence. these must be delimited by '(' and ')', respectively
//leaves the original string unchanged.
std::string remove_ptms(const std::string& sequence);
//converts a sequence into an isosequence. this means replacing all 'I' and 'L' chars with 'J'
//leaves the original string unchanged.
std::string make_iso_sequence(std::string sequence);
//splits a string along the given delimiters. the resulting substrings wil be assigned to the passed vector.
void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ", bool trimEmpty = false);
//a more sophisticated operator== for GenomeCoordinates.
bool same_coordinates(const GenomeCoordinates& lhs, const GenomeCoordinates& rhs);
//compares two genomic coordinates in a way that they can be sorted ascendingly.
//this is important for map/set <GenomeCoordinates> because otherwise they would overwrite 
//other Coordinates even if they arent the same.
bool compare_genome_coordinate_sets_ascending(const std::vector<GenomeCoordinates>& lhs, const std::vector<GenomeCoordinates>& rhs);

//to_string functions for Genomic coordinates.
//as simple string
std::string coordinates_to_string(const GenomeCoordinates& coords);
//as a shortend version.
std::string coordinates_to_short_string(const GenomeCoordinates& coords, unsigned int offset = 0);
//as a line in a gtffile
std::string coordinates_to_gtf_string(const GenomeCoordinates& coords, const std::string& type, bool frameinclude, const std::string& source);
//as a line in a bed file
std::string coordinates_to_bed_string(const GenomeCoordinates& coords, const std::string& name, unsigned int score = 1000);
//as a short bed string
std::string coordinates_to_short_bed_string(const GenomeCoordinates& coords, const std::string& name, unsigned int score = 1000);
//as a line in a gct file
std::string coordinates_to_gct_string(std::vector<GenomeCoordinates> const& coords);
//given a tokenized string this function will generate the resulting genomic coordinates.
GenomeCoordinates extract_coordinates_from_gtf_line(std::vector<std::string> const& tokens);
//this compare_function returns true if lhs.start < rhs.start and false otherwise.
//it should ONLY be used in sort funtions, where only the rough order matters
//for comparisions if two GenomicCoordinates are equal, use compare_genome_coordinate_sets_ascending or compare_coordinates_ascending_whole
bool compare_coordinates_ascending(const GenomeCoordinates& lhs, const GenomeCoordinates& rhs);
//returns true if lhs.start < rhs.start
//otherwise returns true if lhs.end < rhs.end
//otherwise returns false.
bool compare_coordinates_ascending_whole(const GenomeCoordinates& lhs, const GenomeCoordinates& rhs);

//general conveniance to_string method.
template <typename T> std::string to_string(T value) {
	//create an output string stream
	std::ostringstream os;
	//throw the value into the string stream
	os << value;
	//convert the string stream into a string and return
	return os.str();
}

//compares two std::pair<std::string, unsigned int> by their second value.
//returns true if lhs.second is lesser than rhs.second, otherwise returns false.
struct byIntValue {
	bool operator()(const std::pair<std::string, unsigned int>& lhs, const std::pair<std::string, unsigned int>& rhs) const;
};

enum assembly {
	none = 0, 
	primary = 1,
	patchhaploscaff = 2
};

#endif
