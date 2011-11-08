#include <iostream>
#include <iomanip>
#include <boost/bind.hpp>
#include "genfile/SNPDataSink.hpp"
#include "genfile/VCFFormatSNPDataSink.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	namespace {
		std::string int_to_number( std::size_t i ) {
			return string_utils::to_string( i ) ;
		}
	}
	
	VCFFormatSNPDataSink::VCFFormatSNPDataSink( std::string const& filename ):
		m_filename( filename ),
		m_stream_ptr( open_text_file_for_output( filename )),
		m_have_written_header( false ),
		m_number_of_samples( 0 ),
		m_call_threshhold( 0.9 ),
		m_sample_name_getter( &int_to_number )
	{}
	
	void VCFFormatSNPDataSink::write_header( std::size_t number_of_samples ) const {
		(*m_stream_ptr) << "##fileformat=VCFv4.1\n"
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype calls\">\n"
			"##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype call probabilities\">\n"
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			(*m_stream_ptr) << "\t" << m_sample_name_getter( i ) ;
		}
		(*m_stream_ptr) << "\n" ;
	}

	void VCFFormatSNPDataSink::write_snp_impl(
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
		if( !m_have_written_header ) {
			write_header( number_of_samples ) ;
			m_have_written_header = true ;
			m_number_of_samples = number_of_samples ;
		} else if( number_of_samples != m_number_of_samples ) {
			throw genfile::BadArgumentError( "VCFFormatSNPDataSink::write_snp_impl()", "number_of_samples=" + string_utils::to_string( number_of_samples )) ;
		}
		
		char tab = '\t' ;
		(*m_stream_ptr)
			<< chromosome << tab
			<< SNP_position << tab
			<< RSID << tab
			<< first_allele << tab
			<< second_allele << tab
			<< "." << tab
			<< "." << tab
			<< "." << tab
			<< "GT:GP" ;

		std::vector< double > probs( 3 ) ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			probs[0] = get_AA_probability( i ) ;
			probs[1] = get_AB_probability( i ) ;
			probs[2] = get_BB_probability( i ) ;
			std::string call ;
			if( probs[0] >= m_call_threshhold ) {
				call = "0/0" ;
			}
			else if( probs[1] >= m_call_threshhold ) {
				call = "0/1" ;
			}
			else if( probs[2] >= m_call_threshhold ) {
				call = "1/1" ;
			}
			else {
				call = "./." ;
			}
			
			(*m_stream_ptr)
				<< tab << call << ":"
				<< std::fixed << std::setprecision( 6 ) << probs[0] << ","
				<< std::fixed << std::setprecision( 6 ) << probs[1] << ","
				<< std::fixed << std::setprecision( 6 ) << probs[2] ;
		}
		(*m_stream_ptr) << std::endl ;
	}
}

