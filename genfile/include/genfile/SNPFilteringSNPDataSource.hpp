#ifndef GENFILE_SNPFILTERINGSNPDATASOURCE_HPP
#define GENFILE_SNPFILTERINGSNPDATASOURCE_HPP

#include <set>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	class SNPFilteringSNPDataSource: public SNPDataSource
	{
	public:
		// Create a SNPFilteringSNPDataSource from the given source and the given sample indices.
		static std::auto_ptr< SNPFilteringSNPDataSource > create( SNPDataSource& source, std::auto_ptr< SNPIdentifyingDataTest > snp_inclusion_test ) ;
		
	public:
		SNPFilteringSNPDataSource( SNPDataSource& source, std::auto_ptr< SNPIdentifyingDataTest > snp_inclusion_test ) ;

		operator bool() const ;
		unsigned int number_of_samples() const ;
		unsigned int total_number_of_snps() const ;
		unsigned int total_number_of_snps_before_filtering() const ;

		SNPIdentifyingDataTest const& get_snp_inclusion_test() const ;
		
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

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;
		
		void ignore_snp_probability_data_impl() ;



	private:
		
		std::set< std::size_t > get_indices_of_excluded_snps() ;
		
		SNPDataSource& m_source ;
		std::auto_ptr< SNPIdentifyingDataTest > m_snp_inclusion_test ;
		std::set< std::size_t > m_indices_of_excluded_snps ;
	} ;
}

#endif