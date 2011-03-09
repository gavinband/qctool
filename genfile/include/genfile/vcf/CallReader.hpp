#ifndef GENFILE_VCF_CALLREADER_HPP
#define GENFILE_VCF_CALLREADER_HPP

#include <map>
#include <vector>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/function.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/vcf/Types.hpp"

namespace genfile {
	namespace vcf {
		struct CallReader
		//
		// This class reads a set of specific entries of the per-individual data
		// from the given string
		// according.
		// Usage:
		//
		// CallReader( format, data, entry_types )
		//		( "GQ", set_genotype_qualities )
		//		( "GI", set_genotype_intensities ) ;
		//
		// Setters passed as second argument must take a vector of Entry values.
		//
		{
		private:
			typedef boost::function< void ( std::size_t i, std::vector< Entry > const& ) > Setter ;
			typedef std::map< std::size_t, Setter > Setters ;
		private:
			CallReader(
				std::string const& format,
				std::string const& data,
				boost::ptr_map< std::string, VCFEntryType > const& entry_types
			) ;
	
			// Return a vector of pointers to VCFEntryType objects in the supplied map,
			// with the ith VCFEntryType appropriate for the ith format element.
			std::vector< VCFEntryType const* > get_entries_by_position(
				std::vector< std::string > const& format_elts,
				boost::ptr_map< std::string, VCFEntryType > const& entry_types
			) const ;
	
			CallReader& operator()( std::string const& spec, Setter setter ) ;
	
			~CallReader() ;
	
		private:
			std::vector< std::string > const m_format_elts ;
			std::string const& m_data ;
			boost::ptr_map< std::string, VCFEntryType > const& m_entry_types ;
			std::vector< VCFEntryType const* > m_entries_by_position ;
			Setters m_setters ;
			
		private:
			void set_values( std::vector< std::string > const& elts, Setters const& setters ) const ;
			void set_values( std::size_t individual_i, std::string const& elt, Setters const& setters ) const ;
			
			// Forbid copying and assignment.
			CallReader( CallReader const& other ) ;
			CallReader& operator=( CallReader const& other ) ;
		} ;
	}
}

#endif

