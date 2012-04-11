
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNP_DATA_SOURCE_RANDOM_ACCESS_CACHE_HPP
#define SNP_DATA_SOURCE_RANDOM_ACCESS_CACHE_HPP

#include "genfile/SNPDataSource.hpp"
#include "genfile/RandomAccessSNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include <boost/function.hpp>

namespace genfile {
	class SNPDataSourceRandomAccessCache: public RandomAccessSNPDataSource
		// This class reads all of the data from a given SNPDataSource and caches it in memory.
		// Since the quantity of data may be large, this class employs strategies to minimise
		// resident memory usage.
	{
	public:
		typedef std::auto_ptr< SNPDataSourceRandomAccessCache > UniquePtr ;
		typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;

		static UniquePtr create(
			SNPDataSource& source,
			ProgressCallback progress_callback = ProgressCallback()
		) ;
		
		// Construct taking data SNPData
		SNPDataSourceRandomAccessCache(
			SNPDataSource& source,
			ProgressCallback progress_callback = ProgressCallback()
		) ;
		
	public:
		std::size_t get_estimated_memory_usage_in_bytes() const ;
		
		unsigned int number_of_samples() const ; 
		unsigned int number_of_snps() const ;

		// Function: get_snp_identifying_data()
		// Get the SNP ID, RS ID, position, and alleles of the specified snp.
		void get_snp_identifying_data(
			std::size_t snp_index,
			SNPIdentifyingData& snp
		) const ;
		
		// Function: read_snp_probability_data()
		// Read the probability data for the given snp, storing it
		// using the given setter object / function pointer.
		void get_snp_probability_data(
			std::size_t snp_index,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) const ;
		
	private:
		
		void setup( SNPDataSource& source, ProgressCallback progress_callback ) ;

		std::size_t m_number_of_samples ;
		std::vector< SNPIdentifyingData > m_snps ;
		std::vector< std::vector< char > > m_data ;
		mutable std::vector< double > m_genotype_buffer ; // avoid allocating buffer memory on each operation.
		mutable std::vector< char > m_compression_buffer ;
	} ;
}


#endif
