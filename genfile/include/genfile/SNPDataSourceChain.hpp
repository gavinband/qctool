#ifndef SNPDATAsourceCHAIN_HPP
#define SNPDATAsourceCHAIN_HPP

#include <iostream>
#include <string>
#include "../config.hpp"
#if HAVE_BOOST_FUNCTION
#include <boost/function.hpp>
#endif
#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "wildcard.hpp"

namespace genfile {
	// class SNPDataSourceChain represnets a SNPDataSource
	// which gets it data sequentially from a collection of other SNPDataSources
	class SNPDataSourceChain: public SNPDataSource
	{
	public:
		typedef std::auto_ptr< SNPDataSourceChain > UniquePtr ;
		
		// Create a chain of SNPDataSources taking data from the specified files.
		static UniquePtr create(
			std::vector< wildcard::FilenameMatch > const& filenames,
			NotifyProgress notify_progress = NotifyProgress()
		) ;

		// Create a chain of SNPDataSourceRacks taking data from the specified files.
		// Each entry of the outer vector is understood as a list of files which must be chained together.
		// Each entry of the outer vector must therefore have the same size.
		static UniquePtr create(
			std::vector< std::vector< wildcard::FilenameMatch > > const& filenames,
			NotifyProgress notify_progress = NotifyProgress()
		) ;
		
		SNPDataSourceChain() ;
		~SNPDataSourceChain() ;

		void add_source( std::auto_ptr< SNPDataSource > source ) ;
		unsigned int number_of_samples() const ;
		unsigned int total_number_of_snps() const ;
		unsigned int number_of_sources() const ;
		unsigned int number_of_snps_in_source( std::size_t source_index ) const ;
		SNPDataSource const& get_source( std::size_t source_index ) const ;
		operator bool() const ;
		std::string get_source_spec() const ;
		std::string get_summary( std::string const& prefix, std::size_t fill ) const ;

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

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;
		
		void ignore_snp_probability_data_impl() ;

	public:

	#if HAVE_BOOST_FUNCTION
		typedef boost::function< void( int index ) > moved_to_next_source_callback_t ;
	#else
		typedef void( *moved_to_next_source_callback_t )( std::size_t ) ;
	#endif

		void set_moved_to_next_source_callback( moved_to_next_source_callback_t callback ) ;

	private:

		void move_to_next_source() ;
		void move_to_next_nonempty_source_if_necessary() ;

		std::vector< SNPDataSource* > m_sources ;
		std::size_t m_current_source ;
		unsigned int m_number_of_samples ;
		
		moved_to_next_source_callback_t m_moved_to_next_source_callback ;
	} ;
}

#endif