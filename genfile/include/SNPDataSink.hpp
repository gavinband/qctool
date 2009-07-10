#ifndef SNPDATASINK_HPP
#define SNPDATASINK_HPP

#include <iostream>
#include <string>
#include <stddef.h>
#include "snp_data_utils.hpp"
#include "gen.hpp"
#include "bgen.hpp"


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

		virtual operator bool() const { return stream() ; }
		virtual std::ostream& stream() = 0 ;
		virtual std::ostream const & stream() const = 0 ;
		virtual FormatType format() const = 0;

		// Function write_snp().
		// Ideally this would also be a virtual member function.
		// However, a template member function can't be virtual.
		// Therefore, we dispatch to the correct implementation using the format()
		// and stream() members which implementations must provide.
		template<
			typename GenotypeProbabilityGetter
		>
		SNPDataSink& write_snp(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability
		) {
			if( m_number_of_samples == 0 ) {
				m_number_of_samples = number_of_samples ;
			}
			assert(  m_number_of_samples == number_of_samples ) ;

			pre_write_snp() ;

			std::size_t id_field_size = std::min( std::max( SNPID.size(), RSID.size() ), static_cast< std::size_t >( 255 )) ;

			switch( format() ) {
				case e_GenFormat:
					gen::write_snp_block( stream(), number_of_samples, SNPID, RSID, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;
					break ;
				case e_BGenFormat:
					bgen::write_snp_block( stream(), number_of_samples, id_field_size, SNPID, RSID, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;
					break ;
				case e_BGenCompressedFormat:
					bgen::write_compressed_snp_block( stream(), number_of_samples, id_field_size, SNPID, RSID, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;
					break ;
				default:
					assert(0) ; // invalid format type.
			}

			if( stream()) {
				++m_number_of_snps_written ;
				post_write_snp() ;
			}
			return *this ;
		};
	
	protected:
		virtual void pre_write_snp() {} ;
		virtual void post_write_snp() {} ;
		
		uint32_t number_of_samples() const { return m_number_of_samples ; }
		std::size_t number_of_snps_written() const { return m_number_of_snps_written ; }

	private:
	
		uint32_t m_number_of_samples ;
		std::size_t m_number_of_snps_written ;

		SNPDataSink( SNPDataSink const& other ) ;
		SNPDataSink& operator=( SNPDataSink const& other ) ;
	} ;
}

#endif