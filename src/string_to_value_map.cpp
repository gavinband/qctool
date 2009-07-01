#include <string>
#include "string_to_value_map.hpp"

template<>
double string_to_value_map::get_value< double >( std::string const& name ) const {
	return get_double_value( name ) ;
}

template<>
std::string string_to_value_map::get_value< std::string >( std::string const& name ) const {
	return get_string_value( name ) ;
}
