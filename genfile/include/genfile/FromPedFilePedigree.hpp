#ifndef GENFILE_FROM_PED_FILE_PEDIGREE_HPP
#define GENFILE_FROM_PED_FILE_PEDIGREE_HPP

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "genfile/Pedigree.hpp"

namespace genfile {
	class FromPedFilePedigree: public Pedigree
	{
	public:
		FromPedFilePedigree( std::string const& filename ) ;
		FromPedFilePedigree( std::istream& stream ) ;

		std::size_t get_number_of_individuals() const ;
		std::string get_id_of( std::size_t i ) const ;
		
		std::string get_family_of( std::string const& individual ) const ;
		std::vector< std::string > get_parents_of( std::string const& individual ) const ;
		std::vector< std::string > get_children_of( std::string const& individual ) const ;
		Sex get_sex_of( std::string const& individual ) const ;
		
	private:
		std::string m_filename ;
		std::vector< std::string > m_individuals ;
		std::map< std::string, std::string > m_families ;
		std::map< std::string, Sex > m_sexes ;
		std::map< std::string, std::vector< std::string > > m_parents ;
		std::map< std::string, std::vector< std::string > > m_children ;
	private:
		
		void load_data_from_stream( std::istream& stream ) ;
		
		static std::vector< std::string > get_map_value(
			std::map< std::string, std::vector< std::string > > const& map,
			std::string const& key
		) ;
		static std::string get_map_value(
			std::map< std::string, std::string > const& map,
			std::string const& key
		) ;
		static Sex get_map_value(
			std::map< std::string, Sex > const& map,
			std::string const& key
		) ;
	} ;
}


#endif
