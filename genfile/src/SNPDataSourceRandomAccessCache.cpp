
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/SNPDataSource.hpp"
#include "genfile/RandomAccessSNPDataSource.hpp"
#include "genfile/SNPDataSourceRandomAccessCache.hpp"
#include "genfile/zlib.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	namespace impl {
		void compress( std::vector< double > const& genotypes, std::vector< char >* dest ) {
			std::vector< uint16_t > bitbuffer( genotypes.size() ) ;
			for( std::size_t i = 0; i < genotypes.size(); ++i ) {
				double value = std::min( std::max( genotypes[i], 0.0 ), 1.0 ) * 65535.0 ;
				bitbuffer[i] = static_cast< int16_t >( value + 0.5 ) ; // round to nearest int.
			}
			zlib_compress( bitbuffer, dest ) ;
		}

		void uncompress( std::vector< char > const& compressed_data, std::vector< double >* dest ) {
			std::vector< uint16_t > bitbuffer( dest->size() ) ;
			zlib_uncompress( compressed_data, &bitbuffer ) ;
			for( std::size_t i = 0; i < dest->size(); ++i ) {
				(*dest)[i] = static_cast< double >( bitbuffer[i] ) / 65535.0 ;
			}
		}
	}
	
	typedef std::auto_ptr< SNPDataSourceRandomAccessCache > UniquePtr ;

	SNPDataSourceRandomAccessCache::UniquePtr SNPDataSourceRandomAccessCache::create(
		SNPDataSource& source,
		ProgressCallback progress_callback
	) {
		return SNPDataSourceRandomAccessCache::UniquePtr( new SNPDataSourceRandomAccessCache( source, progress_callback )) ;
	}
	
	// Construct taking data SNPData
	SNPDataSourceRandomAccessCache::SNPDataSourceRandomAccessCache(
		SNPDataSource& source,
		ProgressCallback progress_callback
	):
	 	m_number_of_samples( source.number_of_samples() )
	{
		setup( source, progress_callback ) ;
	}
	
	void SNPDataSourceRandomAccessCache::setup(
		SNPDataSource& source,
		ProgressCallback progress_callback
	) {
		m_number_of_samples = source.number_of_samples() ;
		m_genotype_buffer.resize( m_number_of_samples * 3 ) ;

		SNPIdentifyingData snp ;
		while( source.get_snp_identifying_data( snp )) {
			source.read_snp_probability_data( set_genotypes( m_genotype_buffer ) ) ;
			m_snps.push_back( snp ) ;
			impl::compress( m_genotype_buffer, &m_compression_buffer ) ;
			m_data.push_back( m_compression_buffer ) ;
			if( progress_callback ) {
				progress_callback( source.number_of_snps_read(), source.total_number_of_snps() ) ;
			}
		}
		{
			std::vector< char > empty ;
			m_compression_buffer.swap( empty ) ;
		}
		assert( m_data.size() == m_snps.size() ) ;
	}
	
	std::size_t SNPDataSourceRandomAccessCache::get_estimated_memory_usage_in_bytes() const {
		std::size_t result = 0 ;
		for( std::size_t i = 0; i < m_snps.size(); ++i ) {
			result += sizeof( std::size_t ) + sizeof( GenomePosition ) + 2 + m_snps[i].get_SNPID().size() + m_snps[i].get_rsid().size() + m_data[i].capacity() ;
		}
		return result ;
	}

	unsigned int SNPDataSourceRandomAccessCache::number_of_samples() const {
		return m_number_of_samples ;
	}
	
	unsigned int SNPDataSourceRandomAccessCache::number_of_snps() const {
		return m_snps.size() ;
	}

	// Function: get_snp_identifying_data()
	// Get the SNP ID, RS ID, position, and alleles of the specified snp.
	void SNPDataSourceRandomAccessCache::get_snp_identifying_data(
		std::size_t snp_index,
		SNPIdentifyingData& snp
	) const {
		assert( snp_index < m_data.size() ) ;
		snp = m_snps[ snp_index ] ;
	}
	
	// Function: read_snp_probability_data()
	// Read the probability data for the given snp, storing it
	// using the given setter object / function pointer.
	void SNPDataSourceRandomAccessCache::get_snp_probability_data(
		std::size_t snp_index,
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) const {
		assert( snp_index < m_data.size() ) ;
		impl::uncompress( m_data[ snp_index ], &m_genotype_buffer ) ;
		for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
			set_genotype_probabilities(
				i,
				m_genotype_buffer[ (3*i) + 0 ],
				m_genotype_buffer[ (3*i) + 1 ],
				m_genotype_buffer[ (3*i) + 2 ]
			) ;
		}
	}
}
