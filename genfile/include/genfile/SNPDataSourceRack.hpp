#ifndef SNPDATASOURCERACK_HPP
#define SNPDATASOURCERACK_HPP

#include <vector>
#include <memory>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/GenomePosition.hpp"

namespace genfile {
	
	struct SNPDataSourceRack: public SNPDataSource
	{
	public:
		struct Error: public SNPDataError
		{
			Error( std::size_t source_index, SNPIdentifyingData snp )
				: m_source_index( source_index ),
				  m_snp( snp )
			{}

			~Error() throw() {}

			std::size_t source_index() const { return m_source_index ; }
			SNPIdentifyingData snp() const { return m_snp ; }

		private:
			std::size_t m_source_index ;
			SNPIdentifyingData m_snp ;
		} ;

		struct MissingSNPError: public Error
		{
			MissingSNPError( std::size_t source_index, SNPIdentifyingData snp )
				: Error( source_index, snp )
			{}
			// ~MissingSNPError() throw() {}
			char const* what() const throw() { return "MissingSNPError" ; }
			SNPIdentifyingData const& snp() const { return m_snp ; }
		private:
			SNPIdentifyingData const m_snp ;
		} ;

		struct SNPMismatchError: public Error
		{
			SNPMismatchError( std::size_t source_index, SNPIdentifyingData snp )
				: Error( source_index, snp )
			{}
			// ~SNPMismatchError() throw() {}
			char const* what() const throw() { return "SNPMismatchError" ; }
		} ;
		
	public:
		typedef std::auto_ptr< SNPDataSourceRack > UniquePtr ;
		static UniquePtr create( std::vector< wildcard::FilenameMatch > const& filenames ) ;

	public:
		SNPDataSourceRack() ;
		~SNPDataSourceRack() ;
		void add_source( std::auto_ptr< SNPDataSource > source ) ;
		void add_source(
			std::auto_ptr< SNPDataSource > source,
			std::vector< SNPIdentifyingData > const& snps
		) ;
		SNPDataSource& get_source( std::size_t ) const ;

		unsigned int number_of_samples() const ;
		unsigned int total_number_of_snps() const ;
		operator bool() const ;
		std::string get_source_spec() const ;
		std::string get_summary( std::string const& prefix, std::size_t width ) const ;
		std::vector< SNPIdentifyingData > const& get_snps() const { return m_included_snps ; }

	protected:

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
		
		void move_source_to_snp_matching(
			std::size_t source_i,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			char allele1,
			char allele2
		) ;
		
		std::vector< SNPIdentifyingData > get_intersected_snps(
			std::vector< SNPIdentifyingData > const& snps1,
			std::vector< SNPIdentifyingData > const& snps2
		) const ;

		struct RackGenotypeProbabilitySetter
		{
			RackGenotypeProbabilitySetter( GenotypeProbabilitySetter const& base_setter, uint32_t index_of_first_sample ) ;

			void operator()( std::size_t, double, double, double ) const ;
			
		private:
			GenotypeProbabilitySetter const& m_base_setter ;
			uint32_t const m_index_of_first_sample ;
		} ;
		
		std::vector< SNPDataSource* > m_sources ;
		uint32_t m_number_of_samples ;
		std::vector< SNPIdentifyingData > m_included_snps ;
		bool m_read_past_end ;
	} ;
}

#endif
