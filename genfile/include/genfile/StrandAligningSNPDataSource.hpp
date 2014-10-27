
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_STRAND_ALIGNING_SNP_DATA_SOURCE_HPP
#define GENFILE_STRAND_ALIGNING_SNP_DATA_SOURCE_HPP

#include <vector>
#include <map>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	struct StrandAligningSNPDataSource: public SNPDataSource
	{
	public:
		static char const eUnknownStrand = '?' ;
		static char const eForwardStrand = '+' ;
		static char const eReverseStrand  = '-' ;
		static char const eUnknownFlip = '?' ;
		static char const eNoFlip = '+' ;
		static char const eFlip  = '-' ;

		static std::string apply_strand( std::string const& allele, char strand ) ;

	    struct StrandFlipSpec {
	        StrandFlipSpec( char strand_, char flip_ ):
				strand( strand_ ),
	            flip( flip_ )
	        {
	        	assert( strand == eUnknownStrand || strand == eForwardStrand || strand == eReverseStrand ) ;
	        	assert( flip == eUnknownFlip || strand == eFlip || strand == eNoFlip ) ;
	        }

	        StrandFlipSpec():
	            strand( genfile::StrandAligningSNPDataSource::eUnknownStrand ),
	            flip( genfile::StrandAligningSNPDataSource::eUnknownFlip )
	        {}

	        StrandFlipSpec( StrandFlipSpec const& other ):
	            strand( other.strand ),
	            flip( other.flip )
	        {}

	        StrandFlipSpec& operator=( StrandFlipSpec const& other ) {
	            strand = other.strand ;
	            flip = other.flip ;
	            return *this ;
	        }

	        char strand ;
	        char flip ;
	    } ;

	public:
		typedef std::auto_ptr< StrandAligningSNPDataSource > UniquePtr ;
		typedef std::map< genfile::SNPIdentifyingData, StrandFlipSpec, genfile::SNPIdentifyingData::CompareFields > StrandAlignments ; 
		static UniquePtr create( SNPDataSource::UniquePtr source, StrandAlignments const& strand_alignments ) ;

		StrandAligningSNPDataSource( SNPDataSource::UniquePtr source, StrandAlignments const& strand_alignments ) ;

		std::vector< SNPIdentifyingData > const& get_aligned_snps() const { return m_aligned_snps ; }
	private:
		
		SNPDataSource::UniquePtr m_source ;
		StrandAlignments const& m_strand_alignments ;
		std::vector< SNPIdentifyingData > const m_aligned_snps ;
		StrandFlipSpec m_current_strand_flip_spec ;
		
	public:
		operator bool() const { return *m_source ; }
		unsigned int number_of_samples() const { return m_source->number_of_samples() ; }
		OptionalSnpCount total_number_of_snps() const { return m_source->total_number_of_snps() ; }
		std::string get_source_spec() const { return "strand-aligned:" + m_source->get_source_spec() ; }
		SNPDataSource const& get_parent_source() const {
			return *m_source ;
		}
		SNPDataSource const& get_base_source() const {
			return m_source->get_base_source() ;
		}

		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;

	private:

		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;
		StrandFlipSpec get_strand_alignment( SNPIdentifyingData const& snp ) const ;
		VariantDataReader::UniquePtr read_variant_data_impl() ;
		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;
	} ;
}

#endif
