#ifndef SNPDATASINK_HPP
#define SNPDATASINK_HPP

#include <iostream>
#include <string>
#include <stddef.h>
#include <boost/function.hpp>
#include "snp_data_utils.hpp"

namespace genfile {
	struct SNPDataSinkError: public std::exception { char const* what() const throw() { return "SNPDataSinkError" ; } } ;

	// Base class for classes which can store SNP assay data, one snp at a time.
	// After the class is constructed, the intention is that
	// 2. The stream accessible through stream() is writeable so that the snps may be written,
	// one at a time.
	class SNPDataSink: public SNPDataBase
	{
	public:
		SNPDataSink(): m_number_of_samples(0u), m_number_of_snps_written(0u) {} ;
		virtual ~SNPDataSink() {} ;

		// Factory functions
		static std::auto_ptr< SNPDataSink > create( std::string const& filename, std::string const& free_data = "" ) ;
		static std::auto_ptr< SNPDataSink > create( std::string const& filename, CompressionType compression_type, std::string const& free_data = "" ) ;

	public:		

		typedef boost::function< double ( std::size_t ) > GenotypeProbabilityGetter ;

		SNPDataSink& write_snp(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) {
			if( m_number_of_samples == 0 ) {
				m_number_of_samples = number_of_samples ;
			}
			else {
				assert( number_of_samples == m_number_of_samples ) ;
			}
			write_snp_impl( number_of_samples, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;
			if( *this ) {
				++m_number_of_snps_written ;
			}
			return *this ;
		} ;

	public:

		// return the number of samples represented in SNPs in the file.
		// The value returned is undefined until after the first snp has been written.
		uint32_t number_of_samples() const { return m_number_of_samples ; }
		// return the number of SNPs that have been written to the file so far.
		std::size_t number_of_snps_written() const { return m_number_of_snps_written ; }


	public:
		// The following functions must be implemented by derived classes.
		virtual operator bool() const = 0 ;

	protected:
		
		// This function implements the SNP writing, and must be implemented by derived classes.
		virtual void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) = 0 ;

	private:
	
		uint32_t m_number_of_samples ;
		std::size_t m_number_of_snps_written ;

		SNPDataSink( SNPDataSink const& other ) ;
		SNPDataSink& operator=( SNPDataSink const& other ) ;
	} ;
}

#endif