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

		SNPDataSource(): m_number_of_snps_read(0) {} ;
		virtual ~SNPDataSource() {} ;

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

		virtual void read_snp_impl(
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) = 0 ;

		typedef boost::function< bool( std::string const&, std::string const&, uint32_t, char, char ) > SNPMatcher ;

		virtual void read_next_matching_snp_impl(
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2,
			GenotypeProbabilitySetter const& set_genotype_probabilities,
			SNPMatcher const& snp_matcher
		) {
			assert(0) ;
		}

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
		)
		{
			read_snp_impl(
				set_number_of_samples,
				set_SNPID,
				set_RSID,
				set_SNP_position,
				set_allele1,
				set_allele2,
				set_genotype_probabilities
			) ;
			if( *this ) {
				++m_number_of_snps_read ;
			}
			return *this ;
		}
		
		SNPDataSource& read_next_matching_snp(
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities,
			SNPMatcher snp_matcher
		) {
			read_next_matching_snp_impl(
				set_number_of_samples,
				set_SNPID,
				set_RSID,
				set_SNP_position,
				set_allele1,
				set_allele2,
				set_genotype_probabilities,
				snp_matcher
			) ;
			return *this ;
		}
		
		
	private:
		std::size_t m_number_of_snps_read ;
		SNPDataSource( SNPDataSource const& other ) ;
		SNPDataSource& operator=( SNPDataSource const& other ) ;
	} ;
}

#endif