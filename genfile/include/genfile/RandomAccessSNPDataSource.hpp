#ifndef RANDOM_ACCESS_SNP_DATA_SOURCE_HPP
#define RANDOM_ACCESS_SNP_DATA_SOURCE_HPP

#include <boost/function.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	class RandomAccessSNPDataSource
		// This is a base class for classes which represent SNP data (from GEN-like files)
		// which can be accessed by index.
	{
	public:
		virtual ~RandomAccessSNPDataSource() {}

	public:
		typedef std::auto_ptr< RandomAccessSNPDataSource > UniquePtr ;

	public:
		virtual unsigned int number_of_samples() const = 0 ; 
		virtual unsigned int number_of_snps() const = 0 ;

		// Function: get_snp_identifying_data()
		// Get the SNP ID, RS ID, position, and alleles of the specified snp.
		virtual void get_snp_identifying_data(
			std::size_t snp_index,
			SNPIdentifyingData& snp
		) const ;
		
		// Function: read_snp_probability_data()
		// Read the probability data for the given snp, storing it
		// using the given setter object / function pointer.
		virtual void get_snp_probability_data(
			std::size_t snp_index,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) const = 0 ;
	} ;
}

#endif

