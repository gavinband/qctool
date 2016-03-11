
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPDATASINK_HPP
#define SNPDATASINK_HPP

#include <iostream>
#include <string>
#include <map>
#include <stdint.h>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <Eigen/Core>

#include "genfile/snp_data_utils.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/Error.hpp"

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
		typedef boost::shared_ptr< SNPDataSink > SharedPtr ;
		typedef std::multimap< std::string, std::map< std::string, std::string > > Metadata ;
		typedef std::map< std::string, std::vector< VariantEntry > > Info ;
		static std::vector< std::string> get_file_types() ;
		
	public:
		SNPDataSink() ;
		virtual ~SNPDataSink() ;

		// Factory functions
		static UniquePtr create(
			std::string const& filename,
			Metadata const& metadata = Metadata(),
			std::string const& filetype_hint = "guess"
		) ;
	private:
		static UniquePtr create_impl(
			std::string const& filename,
			CompressionType compression_type,
			Metadata const& metadata = Metadata(),
			std::string const& filetype_hint = "guess"
		) ;

	public:		

		typedef boost::function< double ( std::size_t ) > GenotypeProbabilityGetter ;

		typedef boost::function< VariantEntry ( std::size_t ) > SampleNameGetter ;
		SNPDataSink& set_sample_names( std::size_t number_of_samples, SampleNameGetter ) ;
		SNPDataSink& set_metadata( Metadata const& ) ;

		SNPDataSink& write_snp(
			uint32_t number_of_samples,
			VariantIdentifyingData const& snp,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability,
			Info const& info = Info()
		) ;

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
			GenotypeProbabilityGetter const& get_BB_probability,
			Info const& info = Info()
		) ;

		SNPDataSink& write_variant_data(
			VariantIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info = Info()
		) ;
		
		SNPDataSink& finalise() ;

	public:
		// return the number of samples represented in SNPs in the file.
		// The value returned is undefined until after the first snp has been written.
		uint32_t number_of_samples() const { assert( m_samples_have_been_set ) ; return m_number_of_samples ; }
		// return the number of SNPs that have been written to the file so far.
		std::size_t number_of_snps_written() const { return m_number_of_snps_written ; }

		typedef std::pair< SNPDataSink const*, std::ostream::streampos > SinkPos ;
		virtual SinkPos get_stream_pos() const {
			throw OperationUnsupportedError(
				"genfile::SNPDataSink::get_stream_pos()",
				"Get write position",
				"(unknown)"
			) ;
		}
		
		virtual std::string get_spec() const = 0 ;

	public:
		// The following functions must be implemented by derived classes.
		virtual operator bool() const = 0 ;

	protected:
		virtual void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) {} ;
		virtual void set_metadata_impl( Metadata const& ) {} ;
		
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
			GenotypeProbabilityGetter const& get_BB_probability,
			Info const& info
		) ;

		virtual void write_variant_data_impl(
			VariantIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info
		) ;

		virtual void finalise_impl() {}
		
	private:

		uint32_t m_number_of_samples ;
		bool m_samples_have_been_set ;
		std::size_t m_number_of_snps_written ;
		
		Eigen::MatrixXd m_genotypes ;

		SNPDataSink( SNPDataSink const& other ) ;
		SNPDataSink& operator=( SNPDataSink const& other ) ;
	} ;
}

#endif
