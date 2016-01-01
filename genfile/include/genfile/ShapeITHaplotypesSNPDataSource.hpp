
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SHAPEIT_HAPLOTYPES_SNP_DATA_SOURCE_HPP
#define GENFILE_SHAPEIT_HAPLOTYPES_SNP_DATA_SOURCE_HPP

#include <string>
#include <vector>
#include "genfile/IdentifyingDataCachingSNPDataSource.hpp"
#include "genfile/string_utils/slice.hpp"

namespace genfile {
	// Read genotypes from a ShapeIT-style output file.
	// It is assumed that haplotypes come in pairs, corresponding to haplotypes from the same individual.
	class ShapeITHaplotypesSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
	public:
		ShapeITHaplotypesSNPDataSource( std::string const& haplotypes_filename, Chromosome chromosome ) ;
		ShapeITHaplotypesSNPDataSource( std::string const& haplotypes_filename, Chromosome chromosome, CompressionType compression_type ) ;

		Metadata get_metadata() const ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		OptionalSnpCount total_number_of_snps() const { return OptionalSnpCount() ; }
		
		operator bool() const { return m_good ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

		Chromosome chromosome() const { return m_chromosome ; }
		std::string get_source_spec() const { return m_haplotypes_filename ; }

		std::size_t get_number_of_id_columns() const ;
		
		void set_expected_ploidy( GetPloidy ) ;

	private:
		void reset_to_start_impl() ;
		void read_snp_identifying_data_impl( VariantIdentifyingData* result ) ;
		VariantDataReader::UniquePtr read_variant_data_impl() ;
		void ignore_snp_probability_data_impl() ;

	private:
		std::string const m_haplotypes_filename ;
		CompressionType m_compression_type ;
		unsigned int m_number_of_samples ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		bool m_have_chromosome_column ;
		Chromosome m_chromosome ;
		genfile::VariantIdentifyingData m_current_snp ;
		bool m_good ;
		std::string m_current_line ;
		
		GetPloidy m_get_ploidy ;
		std::map< genfile::Chromosome, std::vector< int > > m_ploidies ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void count_samples() ;
		std::vector< int > const& get_or_compute_ploidies( Chromosome const& chromosome ) ;
	} ;
}

#endif
