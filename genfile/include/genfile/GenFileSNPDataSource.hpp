
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILESNPDATAPROVIDER_HPP
#define GENFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "gen.hpp"
#include "IdentifyingDataCachingSNPDataSource.hpp"
#include "vcf/MetadataParser.hpp"

namespace genfile {
	namespace impl {
		struct GenFileSNPDataReader ;
	}
	// This class represents a SNPDataSource which reads its data
	// from a plain GEN file.
	class GenFileSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
		friend struct impl::GenFileSNPDataReader ;
	public:
		GenFileSNPDataSource( std::auto_ptr< std::istream > stream, Chromosome chromosome ) ;
		GenFileSNPDataSource( std::string const& filename, Chromosome chromosome ) ;

		Metadata get_metadata() const ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		OptionalSnpCount total_number_of_snps() const { return m_total_number_of_snps ; }
		
		operator bool() const { return m_stream_ptr->good() ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

		Chromosome chromosome() const { return m_chromosome ; }
		std::string get_source_spec() const { return m_filename ; }

	private:
		void reset_to_start_impl() ;
		
		void read_snp_identifying_data_impl( VariantIdentifyingData* result ) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;

	private:
		std::string m_filename ;
		CompressionType m_compression_type ;
		unsigned int m_number_of_samples ;
		OptionalSnpCount m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		Chromosome m_chromosome ;
		bool m_have_chromosome_column ;
		std::string m_line ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void read_header_data() ;
	} ;
}

#endif
