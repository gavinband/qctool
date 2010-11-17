#ifndef GENFILE_SAMPLEFILTERINGSNPDATASOURCE_HPP
#define GENFILE_SAMPLEFILTERINGSNPDATASOURCE_HPP

#include <set>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"

namespace genfile {
	class SampleFilteringSNPDataSource: public SNPDataSource
	{
	public:
		typedef std::auto_ptr< SampleFilteringSNPDataSource > UniquePtr ;
	public:
		
		// Create a SampleFilteringSNPDataSource from the given source and the given sample indices.
		static std::auto_ptr< SampleFilteringSNPDataSource > create(
			std::auto_ptr< SNPDataSource > source,
			std::set< std::size_t > const& indices_of_samples_to_filter_out
		) ;
		
		SampleFilteringSNPDataSource( std::auto_ptr< SNPDataSource > source, std::set< std::size_t > const& indices_of_samples_to_filter_out ) ;
		~SampleFilteringSNPDataSource() ;

	public:

		operator bool() const ;
		unsigned int number_of_samples() const ;
		unsigned int total_number_of_snps() const ;
		std::string get_source_spec() const ;
		
		SNPDataSource const& get_parent_source() const ;
		SNPDataSource const& get_base_source() const ;

	public:
		struct SampleIndexOutOfRangeError: public SNPDataError
		{
			SampleIndexOutOfRangeError( std::size_t index, std::size_t number_of_samples ): m_index( index ), m_number_of_samples( number_of_samples ) {}
			~SampleIndexOutOfRangeError() throw() {}
			char const* what() const throw() { return "SampleIndexOutOfRangeError" ; }
			std::size_t index_of_sample() const { return m_index ; }
			std::size_t number_of_samples() const { return m_number_of_samples ; }
		private:
			std::size_t const m_index ;
			std::size_t const m_number_of_samples ;
		} ;

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
		void verify_indices( std::set< std::size_t > const& ) const ;
		void read_source_probability_data() ;
		void return_filtered_genotype_probabilities( GenotypeProbabilitySetter const& ) ;
		void set_unfiltered_genotype_probabilities( std::size_t i, double aa, double ab, double bb ) ;

		std::auto_ptr< SNPDataSource > m_source ;
		std::set< std::size_t > m_indices_of_samples_to_filter_out ;
		std::vector< double > m_genotype_data ;
	} ;
}

#endif