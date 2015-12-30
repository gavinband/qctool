
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef LONGFORMATSNPDATASOURCE_HPP
#define LONGFORMATSNPDATASOURCE_HPP

#include <iostream>
#include <string>
#include <map>
#include <boost/unordered_map.hpp>
#include "genfile/SNPDataSource.hpp"

namespace genfile {
	namespace impl {
		struct LongFormatVariantDataReader ;
	}
	
	// This class represents a SNPDataSource which reads its data
	// from a plain GEN file.
	class LongFormatSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
	public:
		LongFormatSNPDataSource( std::string const& filename ) ;
		LongFormatSNPDataSource( std::istream& stream ) ;

	public:
		Metadata get_metadata() const ;

		unsigned int number_of_samples() const { return m_samples.size() ; }
		void get_sample_ids( GetSampleIds ) const ;

		OptionalSnpCount total_number_of_snps() const { return m_variant_map.size() ; }
		operator bool() const { return !m_exhausted ; }

		std::string get_source_spec() const { return m_filename ; }

		std::string const& buffer() const { return m_buffer ; }

	private:
		std::string m_filename ;
		std::size_t const m_buffer_size ;
		std::vector< VariantIdentifyingData > m_variants ;
		typedef std::map< VariantIdentifyingData, int > VariantMap ;
		VariantMap m_variant_map ;
		typedef boost::unordered_map< std::string, int > SampleMap ;
		SampleMap m_sample_map ;
		std::vector< std::string > m_samples ;
		std::string m_buffer ;
		typedef boost::unordered_map< std::pair< int, int >, std::pair< std::size_t, std::size_t > > BufferMap ;
		BufferMap m_buffer_map ;

		std::size_t m_snp_index ;
		bool m_exhausted ;

	private:
		void reset_to_start_impl() ;
		void read_snp_identifying_data_impl( VariantIdentifyingData* result ) ;
		VariantDataReader::UniquePtr read_variant_data_impl() ;
		void ignore_snp_probability_data_impl() ;
		void setup( std::string const& filename ) ;
		void setup( std::istream& stream ) ;
		
		friend struct impl::LongFormatVariantDataReader ;
	} ;
}

#endif
