
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_HAPMAP_HAPLOTYPES_SNP_DATA_SOURCE_HPP
#define GENFILE_HAPMAP_HAPLOTYPES_SNP_DATA_SOURCE_HPP

#include <string>
#include <vector>
#include "genfile/SNPDataSource.hpp"
#include "genfile/string_utils/slice.hpp"

namespace genfile {
	// Read genotypes from a phased haplotypes file in the format available on the HapMap website
	// http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/
	class HapmapHaplotypesSNPDataSource: public SNPDataSource
	{
	public:
		HapmapHaplotypesSNPDataSource( std::auto_ptr< std::istream > stream, Chromosome chromosome ) ;
		HapmapHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome ) ;
		HapmapHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome, CompressionType compression_type ) ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		OptionalSnpCount total_number_of_snps() const { return m_total_number_of_snps ; }
		
		operator bool() const { return *m_stream_ptr ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

		Chromosome chromosome() const { return m_chromosome ; }
		std::string get_source_spec() const { return m_filename ; }

	private:
		void reset_to_start_impl() ;
		
		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;

	private:
		std::string m_filename ;
		CompressionType m_compression_type ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		Chromosome m_chromosome ;
		std::string m_line ;
		std::vector< string_utils::slice > m_elts ;
		std::string m_first_allele, m_second_allele ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void read_header_data() ;
	} ;
}

#endif