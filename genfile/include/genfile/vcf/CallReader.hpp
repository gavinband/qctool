#ifndef GENFILE_VCF_CALLREADER_HPP
#define GENFILE_VCF_CALLREADER_HPP

#include <map>
#include <vector>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/function.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/string_utils/slice.hpp"

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
		{
		private:
			typedef boost::function< void ( std::size_t i, std::vector< Entry > const& ) > Setter ;
			// Setters holds position in the format string (one past the end if field not present),
			// the setter object, and a pointer to the entry type.
			typedef std::multimap< std::size_t, std::pair< Setter, VCFEntryType const* > > Setters ;
		public:
			CallReader(
				std::size_t number_of_samples,
				std::size_t number_of_alleles,
				std::string const& format,
				std::string const& data,
				boost::ptr_map< std::string, VCFEntryType > const& entry_types
			) ;
		
			CallReader& operator()( std::string const& spec, Setter setter ) ;
		
			~CallReader() ;
		
		private:
			std::size_t const m_number_of_samples ;
			std::size_t const m_number_of_alleles ;
			std::vector< std::string > const m_format_elts ;
			std::string const& m_data ;
			boost::ptr_map< std::string, VCFEntryType > const& m_entry_types ;
			std::vector< VCFEntryType const* > m_entries_by_position ;
			Setters m_setters ;
			std::auto_ptr< vcf::GenotypeCallVCFEntryType > m_genotype_entry_type ;
			
		private:
			void set_values( std::vector< string_utils::slice > const& elts, Setters const& setters ) const ;
			void set_values( std::size_t individual_i, string_utils::slice const& elt, Setters const& setters ) const ;
			
			// Forbid copying and assignment.
			CallReader( CallReader const& other ) ;
			CallReader& operator=( CallReader const& other ) ;
		} ;
	}
}

#endif

