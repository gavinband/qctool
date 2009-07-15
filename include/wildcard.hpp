#ifndef GTOOL_WILDCARD_HPP__
#define GTOOL_WILDCARD_HPP__

#include <iostream>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include <GToolException.hpp>

std::pair< std::vector< std::string >, std::vector< std::string > > find_files_matching_path_with_wildcard( std::string filename_with_wildcards, char wildcard_char = '*') ;

#endif

