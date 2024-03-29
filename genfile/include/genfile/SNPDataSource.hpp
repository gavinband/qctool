
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPDATAPROVIDER_HPP
#define SNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "config/config.hpp"
#include <boost/optional.hpp>
#include <boost/function.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/MetadataParser.hpp"
#include "genfile/CohortIndividualSource.hpp"

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
	class SNPDataSource
	{
	public:
		typedef std::auto_ptr< SNPDataSource > UniquePtr ;
		typedef boost::function< void ( std::size_t, std::size_t ) > NotifyProgress ;
		typedef std::multimap< std::string, std::map< std::string, std::string > > Metadata ;
		typedef boost::function< void() > SourceResetCallback ;

	public:
		
		// The following methods are factory functions
		static std::vector< std::string > get_file_types() ;
		static UniquePtr create(
			std::string const& filename,
			Chromosome chromosome_hint = genfile::Chromosome(),
			boost::optional< vcf::MetadataParser::Metadata > const& = boost::optional< vcf::MetadataParser::Metadata >(),
			std::string const& filetype_hint = "guess"
		) ;
	private:
		static UniquePtr create(
			std::string const& filename,
			Chromosome chromosome_hint,
			CompressionType compression_type,
			boost::optional< vcf::MetadataParser::Metadata > const&,
			std::string const& filetype_hint
		) ;

		static UniquePtr create( std::string const& filename, Chromosome, CompressionType compression_type ) ;

	public:
		SNPDataSource() ;
		virtual ~SNPDataSource() ;

		void set_source_reset_callback( SourceResetCallback ) ;
		void reset_to_start() ;

		typedef boost::function< int ( Chromosome const&, std::size_t ) > GetPloidy ;

		virtual Metadata get_metadata() const = 0 ;

		// The next five functions form the main interface for reading snp data. 
		// These typedefs reflect the signatures which the various setter objects
		// needed by these functions must support.
		typedef boost::function< void ( std::size_t, double, double, double ) > GenotypeProbabilitySetter ;
		typedef boost::function< void ( VariantIdentifyingData const& ) > VariantSetter ;
		typedef boost::function< void( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;

		// Function: list_snps
		// Return (via the setter object) a list of all SNPs in the source.
		virtual void list_snps( VariantSetter, ProgressCallback = ProgressCallback() ) ;
		virtual std::vector< VariantIdentifyingData > list_snps( ProgressCallback = ProgressCallback() ) ;

		// Function: read_snp()
		// Read the data for the next snp from the source (and remove it from the source)
		// Store the data using the given setter objects / function pointers.
		// The returned object evaluates to true if the read was successful, otherwise false.
//		SNPDataSource& read_snp(
//			VariantIdentifyingData* variant,
//			GenotypeProbabilitySetter set_genotype_probabilities
//		) ;

		// Function: read_next_snp_with_position_in_range()
		// Move forwards through the source until either
		// 1. a snp with a chromosome/position in the given range is found;
		// 2. a snp with chromosome/position greater than the upper bound is found; or
		// 3. the source is exhausted
		// In case 1., read the snp data and return true; otherwise return false.
		bool get_next_snp_with_position_in_range(
			VariantIdentifyingData* result,
			Chromosome chromosome_lower_bound,
			Chromosome chromosome_upper_bound,
			uint32_t position_lower_bound,
			uint32_t position_upper_bound
		) ;

		// Function: read_next_snp_with_specified_position()
		// Move forwards through the source until either
		// 1. a snp with the given position is found;
		// 2. a snp with position greater than the given position is found; or
		// 3. the source is exhausted
		// In case 1., read the snp data and return true; otherwise return false.
		bool get_next_snp_with_specified_position(
			VariantIdentifyingData* result,
			Chromosome specified_chromosome,
			uint32_t specified_SNP_position
		) ;
		
		// As above but take a GenomePosition.
		bool get_next_snp_with_specified_position(
			VariantIdentifyingData* result,
			GenomePosition specified_position
		) ;
		
		bool get_next_snp_matching(
			VariantIdentifyingData* matching_snp,
			VariantIdentifyingData const& snp,
			VariantIdentifyingData::CompareFields const& comparer = VariantIdentifyingData::CompareFields()
		) ;

		// Function: get_snp_identifying_data()
		// Get identifying data fields and alleles of the next variant in the source.
		// Repeated calls to this function return the data for the same snp, until a call to
		// read_snp_probability_data() or ignore_snp_probability_data() is made.
		SNPDataSource& get_snp_identifying_data(
			VariantIdentifyingData* snp
		) ;

		// Function: read_snp_probability_data()
		// Read the probability data for the next snp in the source, storing it
		// using the given setter object / function pointer.
		// For each snp, you must call get_snp_identifying_data() at least once before
		// calling this function.
		SNPDataSource& read_snp_probability_data(
			GenotypeProbabilitySetter const& set_genotype_probabilities,
			std::string const& genotype_field = ":genotypes:"
		) ;

		// Function: read_variant_data()
		// Return an object that can be used to read data from the next snp in the source.
		VariantDataReader::UniquePtr read_variant_data() ;

		// Function: ignore_snp_probability_data()
		// Read and discard the probability data for the next snp in the source.
		// For each snp, you must call get_snp_identifying_data() at least once before
		// calling this function.
		SNPDataSource& ignore_snp_probability_data() ;
		void pop() { ignore_snp_probability_data() ; }

	public:
		// Return the number of snps which have been read from the source so far.
		std::size_t number_of_snps_read() const { return m_number_of_snps_read ; }

		// Implicit conversion to bool.  Return true if there have been no errors so far.
		virtual operator bool() const = 0 ;
		// Return the number of samples represented in the snps in this source.
		virtual unsigned int number_of_samples() const = 0;

		typedef boost::function< void( std::size_t, std::string ) > GetSampleIds ;
		virtual bool has_sample_ids() const = 0 ;
		virtual void get_sample_ids( GetSampleIds ) const ;

		// Return the total number of snps the source contains.
		typedef boost::optional< std::size_t > OptionalSnpCount ;
		virtual OptionalSnpCount total_number_of_snps() const = 0 ; // deprecated
		OptionalSnpCount size() const { return total_number_of_snps() ; }

		// Return a string identifying the source of the SNP data
		virtual std::string get_source_spec() const = 0 ;
		
		// In order to support filtered views of source data, we interpret the source
		// as forming a hierarchy.  This method returns the parent in the hierarchy.
		// For the base class, we just return this object.
		virtual SNPDataSource const& get_parent_source() const {
			return *this ;
		}

		// This method returns the base (i.e. eventual parent) in the hierarchy.
		// For the base class, we just return this object.
		virtual SNPDataSource const& get_base_source() const {
			return *this ;
		}
		
		virtual std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	protected:

		virtual void get_snp_identifying_data_impl(
			VariantIdentifyingData* variant
		) = 0 ;

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& ,
			std::string const& genotype_field
		) ;

		virtual VariantDataReader::UniquePtr read_variant_data_impl() = 0 ;

		virtual void ignore_snp_probability_data_impl() = 0 ;
		virtual void reset_to_start_impl() = 0 ;

	protected:

		enum State { e_HaveReadIdentifyingData, e_HaveNotReadIdentifyingData } ;
		State const& state() const { return m_state ; }

	private:
		std::size_t m_number_of_snps_read ;
		SNPDataSource( SNPDataSource const& other ) ;
		SNPDataSource& operator=( SNPDataSource const& other ) ;

		// state variable SNP identifying data
		State m_state ;
		
		SourceResetCallback m_source_reset_callback ;
	} ;
}

#endif
