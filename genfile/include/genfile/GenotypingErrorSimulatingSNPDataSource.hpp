#ifndef GENFILE_GENOTYPINGERRORSIMULATINGSNPDATASOURCE_HPP
#define GENFILE_GENOTYPINGERRORSIMULATINGSNPDATASOURCE_HPP

#include <set>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"

namespace genfile {
	class GenotypingErrorSimulatingSNPDataSource: public SNPDataSource
	{
	public:
		static std::auto_ptr< GenotypingErrorSimulatingSNPDataSource > create(
			std::auto_ptr< SNPDataSource > source,
			std::set< std::size_t > const& indices_of_error_snps
		) ;

		static std::auto_ptr< GenotypingErrorSimulatingSNPDataSource > create(
			std::auto_ptr< SNPDataSource > source,
			double probability_of_error
		) ;
		
		GenotypingErrorSimulatingSNPDataSource( std::auto_ptr< SNPDataSource > source, std::set< std::size_t > const&  positions_of_error_snps ) ;
		~GenotypingErrorSimulatingSNPDataSource() ;

	public:

		unsigned int number_of_samples() const ;
		unsigned int total_number_of_snps() const ;
		operator bool() const ;

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

		void ignore_snp_probability_data_impl() ;

	private:
		void verify_indices( std::set< std::size_t > const& ) const ;
		void read_source_probability_data() ;
		void return_filtered_genotype_probabilities( GenotypeProbabilitySetter const& ) ;
		void set_source_genotype_probabilities( std::size_t i, double aa, double ab, double bb ) ;

		SNPDataSource& m_source ;
		std::vector< double > m_genotype_data ;
	} ;
}

#endif