#ifndef SNPDATASINK_HPP
#define SNPDATASINK_HPP

#include <iostream>
#include <string>
#include <stdint.h>
#include <boost/function.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	struct SNPDataSinkError: public std::exception { char const* what() const throw() { return "SNPDataSinkError" ; } } ;

	// Base class for classes which can store SNP assay data, one snp at a time.
	// After the class is constructed, the intention is that
	// 2. The stream accessible through stream() is writeable so that the snps may be written,
	// one at a time.
	class SNPDataSink
	{
	public:
		typedef std::auto_ptr< SNPDataSink > UniquePtr ;
	public:
		SNPDataSink() ;
		virtual ~SNPDataSink() ;

		// Factory functions
		static UniquePtr create( std::string const& filename, std::string const& free_data = "" ) ;
	private:
		static UniquePtr create_impl( std::string const& filename, CompressionType compression_type, std::string const& free_data = "" ) ;

	public:		

		typedef boost::function< double ( std::size_t ) > GenotypeProbabilityGetter ;

		SNPDataSink& write_snp(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) ;

		SNPDataSink& write_variant_data(
			SNPIdentifyingData const& id_data,
			VariantDataReader& data_reader
		) ;

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
			std::string first_allele,
			std::string second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) = 0 ;

		virtual void write_variant_data_impl(
			SNPIdentifyingData const& id_data,
			VariantDataReader& data_reader
		) ;

	private:
	
		uint32_t m_number_of_samples ;
		std::size_t m_number_of_snps_written ;

		SNPDataSink( SNPDataSink const& other ) ;
		SNPDataSink& operator=( SNPDataSink const& other ) ;
	} ;
}

#endif
