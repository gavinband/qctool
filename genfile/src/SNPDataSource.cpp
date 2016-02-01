
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <boost/optional.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/GenFileSNPDataSource.hpp"
#include "genfile/BGenFileSNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/VCFFormatSNPDataSource.hpp"
#include "genfile/HapmapHaplotypesSNPDataSource.hpp"
#include "genfile/ImputeHaplotypesSNPDataSource.hpp"
#include "genfile/ShapeITHaplotypesSNPDataSource.hpp"
#include "genfile/DosageFileSNPDataSource.hpp"
#include "genfile/BedFileSNPDataSource.hpp"
#include "genfile/LongFormatSNPDataSource.hpp"
#include "genfile/HLAIMPAsBiallelicVariantDataSource.hpp"
#include "genfile/get_set.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/MetadataParser.hpp"

namespace genfile {
	std::vector< std::string > SNPDataSource::get_file_types() {
		std::vector< std::string > result ;
		result.push_back( "gen" ) ;
		result.push_back( "bgen" ) ;
		result.push_back( "bgen_v12" ) ;
		result.push_back( "vcf" ) ;
		result.push_back( "hapmap_haplotypes" ) ;
		result.push_back( "impute_haplotypes" ) ;
		result.push_back( "shapeit_haplotypes" ) ;
		result.push_back( "binary_ped" ) ;
		result.push_back( "long" ) ;
		result.push_back( "hlaimp" ) ;
		return result ;
	}
	
	std::auto_ptr< SNPDataSource > SNPDataSource::create(
		std::string const& filename,
		Chromosome chromosome_hint,
		boost::optional< vcf::MetadataParser::Metadata > const& metadata,
		std::string const& filetype_hint
	) {
		std::pair< std::string, std::string > uf = uniformise( filename ) ;
		
		if( filetype_hint != "guess" ) {
			uf.first = filetype_hint ;
		}
		
		if( uf.first == "bgen" ) {
			return std::auto_ptr< SNPDataSource >( new BGenFileSNPDataSource( uf.second, chromosome_hint )) ;
		}
		else if( uf.first == "vcf" ) {
			return SNPDataSource::UniquePtr( new VCFFormatSNPDataSource( uf.second, metadata )) ;
		}
		else if( uf.first == "gen" ) {
			return std::auto_ptr< SNPDataSource >( new GenFileSNPDataSource( uf.second, chromosome_hint )) ;
		}
		else if( uf.first == "dosage" ) {
			return std::auto_ptr< SNPDataSource >( new DosageFileSNPDataSource( uf.second, chromosome_hint )) ;
		}
		else if( uf.first == "hapmap_haplotypes" ) {
			return std::auto_ptr< SNPDataSource >( new HapmapHaplotypesSNPDataSource( uf.second, chromosome_hint )) ;
		}
		else if( uf.first == "impute_haplotypes" ) {
			return std::auto_ptr< SNPDataSource >( new ImputeHaplotypesSNPDataSource( uf.second, chromosome_hint )) ;
		}
		else if( uf.first == "shapeit_haplotypes" ) {
			return std::auto_ptr< SNPDataSource >( new ShapeITHaplotypesSNPDataSource( uf.second, chromosome_hint )) ;
		}
		else if( uf.first == "binary_ped" ) {
			if( uf.second.size() < 4 || uf.second.substr( uf.second.size() - 4, 4 ) != ".bed" ) {
				throw genfile::BadArgumentError(
					"SNPDataSource::create()",
					"filename=\"" + uf.second + "\"",
					"For binary PED format, expected the .bed extension."
				) ;
			}
			std::string const bedFilename = uf.second ;
			std::string bimFilename = bedFilename ;
			bimFilename.replace( bimFilename.size() - 4, 4, ".bim" ); 
			std::string famFilename = bedFilename ;
			famFilename.replace( famFilename.size() - 4, 4, ".fam" ); 
			return std::auto_ptr< SNPDataSource >( new BedFileSNPDataSource(
				uf.second, bimFilename, famFilename
			) ) ;
		}
		else if( uf.first == "long" ) {
			return std::auto_ptr< SNPDataSource >( new LongFormatSNPDataSource( uf.second ) ) ;
		}
		else if( uf.first == "hlaimp" ) {
			return std::auto_ptr< SNPDataSource >( new HLAIMPAsBiallelicVariantDataSource( uf.second ) ) ;
		}
		else {
			throw genfile::BadArgumentError(
				"genfile::SNPDataSource::create()",
				"filetype_hint=\"" + filetype_hint + "\"",
				"Unrecognised file type."
			) ;
			// assume GEN format.
			// return std::auto_ptr< SNPDataSource >( new GenFileSNPDataSource( uf.second, chromosome_hint, compression_type )) ;
		}
	}

	SNPDataSource::SNPDataSource():
		m_number_of_snps_read(0),
		m_state( e_HaveNotReadIdentifyingData ),
		m_source_reset_callback(0)
	{}

	SNPDataSource::~SNPDataSource() {}

	void SNPDataSource::reset_to_start() {
		reset_to_start_impl() ;
		m_number_of_snps_read = 0 ;
		m_state = e_HaveNotReadIdentifyingData ;
		if( m_source_reset_callback ) {
			m_source_reset_callback() ;
		}
	}

	void SNPDataSource::set_source_reset_callback( SourceResetCallback callback ) {
		m_source_reset_callback = callback ;
	}

	void SNPDataSource::list_snps( VariantSetter setter, boost::function< void( std::size_t, boost::optional< std::size_t > ) > progress_callback ) {
		reset_to_start() ;
		for(
			VariantIdentifyingData variant ;
			get_snp_identifying_data( &variant ) ;
			ignore_snp_probability_data()
		) {
			setter( variant ) ;
			if( progress_callback ) {
				progress_callback( number_of_snps_read() + 1, total_number_of_snps() ) ;
			}
		}
		
	}

	std::vector< VariantIdentifyingData > SNPDataSource::list_snps( ProgressCallback callback ) {
		std::vector< VariantIdentifyingData > result ;
		list_snps( boost::bind( &std::vector< VariantIdentifyingData >::push_back, &result, _1 ), callback ) ;
		return result ;
	}

/*	SNPDataSource& SNPDataSource::read_snp(
		VariantIdentifyingData* result,
		GenotypeProbabilitySetter set_genotype_probabilities
	)
	{
		VariantIdentifyingData variant ;
		get_snp_identifying_data( &variant ) ;
		if( *this ) {
			*result = variant ;
			read_snp_probability_data( set_genotype_probabilities ) ;
		}
		return *this ;
	}
*/
	bool SNPDataSource::get_next_snp_with_specified_position(
		VariantIdentifyingData* result,
		Chromosome specified_chromosome,
		uint32_t specified_SNP_position
	) {
		return get_next_snp_with_position_in_range(
			result,
			specified_chromosome,
			specified_chromosome,
			specified_SNP_position,
			specified_SNP_position
		) ;
	}

	bool SNPDataSource::get_next_snp_with_specified_position(
		VariantIdentifyingData* result,
		GenomePosition specified_position
	) {
		return get_next_snp_with_position_in_range(
			result,
			specified_position.chromosome(),
			specified_position.chromosome(),
			specified_position.position(),
			specified_position.position()
		) ;
	}

	bool SNPDataSource::get_next_snp_with_position_in_range(
		VariantIdentifyingData* result,
		Chromosome chromosome_lower_bound,
		Chromosome chromosome_upper_bound,
		uint32_t position_lower_bound,
		uint32_t position_upper_bound
	) {
		assert( chromosome_lower_bound <= chromosome_upper_bound ) ;
		assert( position_lower_bound <= position_upper_bound ) ;

		VariantIdentifyingData variant ;
		while( *this ) {
			get_snp_identifying_data( &variant ) ;
			if( *this ) {
				genfile::Chromosome const& chromosome = variant.get_position().chromosome() ;
				if( chromosome < chromosome_lower_bound || ((chromosome <= chromosome_upper_bound) && (variant.get_position().position() < position_lower_bound ))) {
					ignore_snp_probability_data() ;
				}
				else if( chromosome > chromosome_upper_bound || variant.get_position().position() > position_upper_bound ) {
					return false ;
				}
				else {
					*result = variant ;
					return *this ;
				}
			}
		}

		return false ;	
	}

	bool SNPDataSource::get_next_snp_matching(
		VariantIdentifyingData* result,
		VariantIdentifyingData const& variant_to_match,
		VariantIdentifyingData::CompareFields const& comparer
	) {
		assert( result ) ;
		bool found = false ;
		VariantIdentifyingData variant ;
		while(
			get_next_snp_with_specified_position(
				&variant,
				variant_to_match.get_position()
			)
		) {
			if( comparer.are_equal( variant, variant_to_match )) {
				found = true ;
				break ;
			} else {
				ignore_snp_probability_data() ;
			}
		}
		if( found ) {
			*result = variant ;
		}
		return found ;
	}

	SNPDataSource& SNPDataSource::get_snp_identifying_data(
		VariantIdentifyingData* result
	) {
		VariantIdentifyingData variant ;
		get_snp_identifying_data_impl( &variant ) ;
		if( *this ) {
			m_state = e_HaveReadIdentifyingData ;
			*result = variant ;
		}
		return *this ;
	}

	SNPDataSource& SNPDataSource::read_snp_probability_data(
		GenotypeProbabilitySetter const& set_genotype_probabilities,
		std::string const& genotype_field
	) {
		assert( m_state == e_HaveReadIdentifyingData ) ;
		read_snp_probability_data_impl( set_genotype_probabilities, genotype_field ) ;
		if( *this ) {
			m_state = e_HaveNotReadIdentifyingData ;
			++m_number_of_snps_read ;
		} else {
			throw MalformedInputError( get_source_spec(), m_number_of_snps_read ) ;
		}
		return *this ;
	}

 	void SNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities,
		std::string const& genotype_field
	) {
		VariantDataReader::UniquePtr reader = read_variant_data_impl() ;
		if( reader.get() ) {
			vcf::GenotypeSetter< GenotypeProbabilitySetter > setter( set_genotype_probabilities ) ;
			reader->get( genotype_field, setter ) ;
		}
	}
	
	VariantDataReader::UniquePtr SNPDataSource::read_variant_data() {
		assert( m_state == e_HaveReadIdentifyingData ) ;
		VariantDataReader::UniquePtr result = read_variant_data_impl() ;
		if( *this ) {
			m_state = e_HaveNotReadIdentifyingData ;
			++m_number_of_snps_read ;
		}
		return result ;
	}
	
	SNPDataSource& SNPDataSource::ignore_snp_probability_data() {
		assert( m_state == e_HaveReadIdentifyingData ) ;
		ignore_snp_probability_data_impl() ;
		if( *this ) {
			m_state = e_HaveNotReadIdentifyingData ;
			++m_number_of_snps_read ;
		}
		return *this ;
	}

	std::string SNPDataSource::get_summary( std::string const& prefix, std::size_t width ) const {
		std::ostringstream ostr ;
		ostr << prefix << std::setw( width ) << "Spec: " << get_source_spec() << "\n" ;
		ostr << prefix << std::setw( width ) << "Number of samples: " <<  number_of_samples() << "\n" ;
		ostr << prefix << std::setw( width ) << "Number of SNPs: " <<  total_number_of_snps() << "\n" ;
		return ostr.str() ;
	}	

}

