#if !defined	MAIN_FUNCTIONS_H
#define			MAIN_FUNCTIONS_H

#include "Parsers\GTFParser.h"

//------------------------------------------------------
//this header was intruduced to declutter the main.cpp file.
//------------------------------------------------------

inline bool cmd_option_exists(char** begin, char** end, const std::string& option) {
	return std::find(begin, end, option) != end;
}


// parameter_is_set checks whether a parameter is set in the command line arguments or not
// If set: The value is saved in the 'value' string which is returned by the reference and the function returns true
// Else:   The value will be ignored and the 'value' stays the same. At the same time the function retuns false
inline bool parameter_is_set(std::vector<std::pair<std::string, std::string>> const& args, std::string const& key, std::string& value) {
	bool param_is_set = false;
	for (size_t j(0); j < args.size(); ++j) {
		if (args[j].first == key) {
			value = args[j].second;
			param_is_set = true;
			break;
		}
	}
	return param_is_set;
}


// This function converts strings to lower case. It is used to convert the command line argument keys to lower case.
// It enables the opportunity to ignore case sensitivity regarding the command line argument keys.
inline void to_lowercase(std::string& str) {
	for (size_t i(0); i < str.size(); i++) {
		str[i] = tolower(str[i]);
	}
}


//Command line argument parser that saves odd m_entries as keys and the following even m_entries as their values.
// in a std::map<std::string, std::string> the command line arguments are finally saves as key-value pairs
inline std::vector<std::pair<std::string, std::string>> get_argument_list(int argc, char* argv[]) {
	std::vector<std::pair<std::string, std::string>> arguments = std::vector<std::pair<std::string, std::string>>();
	std::string key = "";
	for (int i(1); i < argc; ++i) {
		if (i % 2 == 1) { // key which contains usually -command like -q, -lib, -o, ...
			key = argv[i];
			to_lowercase(key); //converts the keys to lower case, i.e.: MSMSTolerance to msmstolerance
		} else { // value like settings or filepaths
			arguments.push_back(std::pair<std::string, std::string>(key, argv[i]));
		}
	}
	return arguments;
}


// The print_parameter_list function prints all the possible command line argument keys and the related information (formatted).
// This way some kind of help menu can be displayed to describe the command line arguments.
inline void print_parameter_list(std::vector<std::pair<std::string, std::string>> paramList) {
	std::cout << "\n";
	for (size_t i(0); i < paramList.size(); ++i) {
		std::string first = paramList[i].first;
		std::string second = paramList[i].second;
		if (first == "-q" || first == "-o") {
			int no = 20;
			std::cout << std::left << std::setw(no) << first << second << std::endl << std::endl;
		} else {
			int no = 20;
			if (second.size() < 60) {
				std::cout << std::left << std::setw(no) << first << second << "\n\n";
			} else {
				std::cout << std::left << std::setw(no) << first;
				while (second.size() > 60) {
					int last_space_pos = 55;
					for (size_t j = 55; j >= 0; j--) {
						if (second[j] == ' ') {
							last_space_pos = j;
							break;
						}
					}
					for (int j = 0; j <= last_space_pos; ++j) {
						std::cout << second[j];
					}
					std::cout << std::endl;
					second = second.substr(last_space_pos + 1, second.size());
					std::cout << std::left << std::setw(no) << " ";
				}
				std::cout << second << std::endl << std::endl;
			}
		}
	}
	std::cout << "\n";
}

#endif
