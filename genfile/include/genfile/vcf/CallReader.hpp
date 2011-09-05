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
#include "genfile/VariantDataReader.hpp"

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
		public:
			// typedef boost::function< void ( std::size_t i, std::vector< Entry >& ) > Setter ;
			typedef VariantDataReader::PerSampleSetter Setter ;
		public:
			typedef std::auto_ptr< CallReader > UniquePtr ;
			
			CallReader(
				std::size_t number_of_samples,
				std::size_t number_of_alleles,
				std::string const& format,
				std::string const& data,
				boost::ptr_map< std::string, VCFEntryType > const& entry_types
			) ;
		
			CallReader& get( std::string const& spec, Setter& setter ) ;
		
			void get_format_elts( boost::function< void ( std::string ) > ) const ;
		
			~CallReader() {} ;
		
		private:
			std::size_t const m_number_of_samples ;
			std::size_t const m_number_of_alleles ;
			std::vector< std::string > const m_format_elts ;
			std::string const m_data ;
			boost::ptr_map< std::string, VCFEntryType > const& m_entry_types ;
			std::size_t m_index ;
			std::vector< VCFEntryType const* > m_entries_by_position ;
			GenotypeCallVCFEntryType m_genotype_call_entry_type ;
			std::vector< string_utils::slice > m_components ;
			std::vector< std::size_t > m_component_counts ;
			// std::vector< std::size_t > m_components_sizes ;
			// std::vector< string_utils::slice > m_components ;
			std::vector< Setter::Integer > m_genotype_calls ;
			std::vector< std::size_t > m_ploidy ;
		private:
			void split_data() ;
			void load_genotypes() ;
			void set_values(
				std::size_t individual_i,
				string_utils::slice const* begin_components,
				string_utils::slice const* end_components,
				std::size_t element_pos,
				VCFEntryType const& entry_type,
				Setter& setter
			) ;

			void unsafe_set_values(
				std::size_t individual_i,
				string_utils::slice const* begin_components,
				string_utils::slice const* end_components,
				std::size_t element_pos,
				VCFEntryType const& entry_type,
				Setter& setter
			) ;
			
			// Forbid copying and assignment.
			CallReader( CallReader const& other ) ;
			CallReader& operator=( CallReader const& other ) ;
		} ;
	}
}

#endif

