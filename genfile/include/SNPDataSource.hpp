#ifndef SNPDATAPROVIDER_HPP
#define SNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include "snp_data_utils.hpp"

namespace genfile {
	struct SNPDataSourceError: public SNPDataError { char const* what() const throw() { return "SNPDataSourceError" ; } } ;
	struct FileIsInvalidError: public SNPDataSourceError { char const* what() const throw() { return "FileIsInvalidError" ; } } ;
	struct FileHasTrailingDataAfterLastSNP: public FileIsInvalidError { char const* what() const throw() { return "FileHasTrailingDataAfterLastSNP" ; } } ;
	struct FileContainsSNPsOfDifferentSizes: public FileIsInvalidError { char const* what() const throw() { return "FileContainsSNPsOfDifferentSizes" ; } } ;

	// Base class for classes which provide SNP assay data, one snp at a time.
	// After the class is constructed, the intention is that
	// 1. the number_of_samples() and total_number_of_snps() functions return information
	// reflecting the data in the file or source, and
	// 2. The stream accessible through stream() is readable and pointing to the first snp in the file.
	class SNPDataSource: public SNPDataBase
	{
	public:

		SNPDataSource() ;
		virtual ~SNPDataSource() ;

		// The following methods must be overriden in derived classes
	public:
		virtual operator bool() const = 0 ;
		virtual unsigned int number_of_samples() const = 0;
		virtual unsigned int total_number_of_snps() const = 0 ;

	public:
		std::size_t number_of_snps_read() const { return m_number_of_snps_read ; }

	protected:

		typedef boost::function< void ( uint32_t ) > IntegerSetter ;
		typedef boost::function< void ( std::string const& ) > StringSetter ;
		typedef boost::function< void ( char ) > AlleleSetter ;
		typedef boost::function< void ( uint32_t ) > SNPPositionSetter ;
		typedef boost::function< void ( std::size_t, double, double, double ) > GenotypeProbabilitySetter ;

		virtual void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) = 0 ;	

		virtual void read_snp_probability_data_impl(
		 	uint32_t* number_of_samples,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) = 0 ;

		virtual void ignore_snp_probability_data_impl(
			uint32_t number_of_samples
		) = 0 ;

		virtual void read_snp_impl(
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;

		typedef boost::function< bool( std::string const&, std::string const&, uint32_t, char, char ) > SNPMatcher ;

		virtual std::size_t read_next_matching_snp_impl(
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2,
			GenotypeProbabilitySetter const& set_genotype_probabilities,
			SNPMatcher const& snp_matcher
		) ;

	public:
		// The following methods are factory functions
		static std::auto_ptr< SNPDataSource > create( std::string const& filename ) ;
		static std::auto_ptr< SNPDataSource > create( std::string const& filename, CompressionType compression_type ) ;
		static std::auto_ptr< SNPDataSource > create( std::vector< std::string > const& filenames ) ;

		// Function read_snp().
		// This is the method which returns snp data from the source.
		// Ideally this would also be a virtual member function.
		// However, a template member function can't be virtual.
		// Therefore, we dispatch to the correct implementation using the format()
		// and stream() members which implementations must provide.
		SNPDataSource& read_snp(
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities
		) ;
		
		SNPDataSource& read_next_matching_snp(
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities,
			SNPMatcher snp_matcher
		) ;
		
		virtual void get_snp_identifying_data(
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;

		virtual void read_snp_probability_data(
			uint32_t* number_of_samples,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;

		virtual void ignore_snp_probability_data( uint32_t number_of_samples ) ;

	protected:

		enum State { e_HaveReadIdentifyingData, e_HaveNotReadIdentifyingData } ;
		State const& state() const { return m_state ; }

	private:
		std::size_t m_number_of_snps_read ;
		SNPDataSource( SNPDataSource const& other ) ;
		SNPDataSource& operator=( SNPDataSource const& other ) ;

		// state variable SNP identifying data
		State m_state ;
	} ;

	class IdentifyingDataCachingSNPDataSource: public SNPDataSource
	{
		virtual void read_snp_identifying_data_impl( 
			uint32_t* number_of_samples,
			std::string* SNPID,
			std::string* RSID,
			uint32_t* SNP_position,
			char* allele1,
			char* allele2
		) = 0 ;

		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;

	private:
		bool m_have_cached_identifying_data ;
		uint32_t m_cached_number_of_samples, m_cached_SNP_position ;
		std::string m_cached_SNPID, m_cached_RSID ;
		char m_cached_allele1, m_cached_allele2 ;
	} ;
}

#endif