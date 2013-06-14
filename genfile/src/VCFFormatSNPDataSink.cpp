
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
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
		m_compression_type( get_compression_type_indicated_by_filename( filename ) ),
		m_stream_ptr( open_text_file_for_output( filename, m_compression_type )),
		m_have_written_header( false ),
		m_number_of_samples( 0 ),
		m_call_threshhold( 0.9 )
	{
		(*m_stream_ptr) << std::resetiosflags( std::ios::floatfield ) << std::setprecision( 6 ) ;
	}
	
	std::string VCFFormatSNPDataSink::get_spec() const {
		return m_filename ;
	}
	
	void VCFFormatSNPDataSink::write_header( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) const {
		(*m_stream_ptr) << "##fileformat=VCFv4.1\n"
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype calls\">\n"
			//"##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype call probabilities\">\n"
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			(*m_stream_ptr) << "\t" << sample_name_getter( i ) ;
		}
		(*m_stream_ptr) << "\n" ;
	}

	namespace {
		struct DataWriter: public VariantDataReader::PerSampleSetter {
			DataWriter( boost::ptr_vector< std::ostringstream >& streams ):
				m_streams( streams )
			{}
			
			~DataWriter() throw() {}

			void set_number_of_samples( std::size_t n ) {
				assert( m_streams.size() == n ) ;
			}

			void set_sample( std::size_t i ) {
				m_sample_i = i ;
				m_streams[i].str( std::string() ) ;
			}

			void set_number_of_entries( std::size_t n ) {
				m_number_of_entries = n ;
			}

			void set_order_type( OrderType const type ) {
				m_order_type == type ;
			}

			void operator()( MissingValue const value ) {
				m_streams[ m_sample_i ] << (( m_entry_i > 0 ) ? "," : "" ) << value ;
			}

			void operator()( std::string& value ) {
				m_streams[ m_sample_i ] << (( m_entry_i > 0 ) ? "," : "" ) << value ;
			}

			void operator()( Integer const value ) {
				m_streams[ m_sample_i ] << (( m_entry_i > 0 ) ? "," : "" ) << value ;
			}

			void operator()( double const value ) {
				m_streams[ m_sample_i ] << (( m_entry_i > 0 ) ? "," : "" ) << value ;
			}

		private:
			boost::ptr_vector< std::ostringstream >& m_streams ;
			std::size_t m_sample_i ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
			OrderType m_order_type ;
		} ;
	}

	void VCFFormatSNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		char const tab = '\t' ;
		(*m_stream_ptr)
			<< id_data.get_position().chromosome() << tab
			<< id_data.get_position().position() << tab
			<< id_data.get_rsid() << tab
			<< id_data.get_first_allele() << tab
			<< id_data.get_second_allele() << tab
			<< "." << tab
			<< "." << tab ;

		write_info( info ) ;

		(*m_stream_ptr) << tab << "GT" ;//:GP" ;

		DataWriter writer( m_streams) ;
	}

/*
	void VCFFormatSNPDataSink::write_snp_impl(
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
	) {
		if( number_of_samples != m_number_of_samples ) {
			throw genfile::BadArgumentError( "VCFFormatSNPDataSink::write_snp_impl()", "number_of_samples=" + string_utils::to_string( number_of_samples )) ;
		}
		char const tab = '\t' ;
		(*m_stream_ptr)
			<< chromosome << tab
			<< SNP_position << tab
			<< RSID << tab
			<< first_allele << tab
			<< second_allele << tab
			<< "." << tab
			<< "." << tab ;

		write_info( info ) ;

		(*m_stream_ptr) << tab << "GT" ;//:GP" ;
		
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
				<< tab << call ;
				//<< ":"
				//<< probs[0] << ","
				//<< probs[1] << ","
				//<< probs[2] ;
				
		}
		(*m_stream_ptr) << std::endl ;
	}
	*/
	void VCFFormatSNPDataSink::write_info( Info const& info ) {
		if( info.empty() ) {
			(*m_stream_ptr) << "." ;
		}
		else {
			Info::const_iterator info_i = info.begin(), end_info_i = info.end() ;
			for( int info_count = 0; info_i != end_info_i; ++info_i, ++info_count ) {
				if( info_count > 0 ) {
					(*m_stream_ptr) << ";" ;
				}
				(*m_stream_ptr) << info_i->first << "=" ;
				for( std::size_t i = 0; i < info_i->second.size(); ++i ) {
					if( i > 0 ) {
						(*m_stream_ptr) << "," ;
					}
					if( info_i->second[i].is_missing() ) {
						(*m_stream_ptr) << "." ;
					}
					else {
						(*m_stream_ptr) << info_i->second[i] ;
					}
				}
			}
		}
	}

	void VCFFormatSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) {
		assert( sample_name_getter ) ;
		assert( !m_have_written_header ) ;
		write_header( number_of_samples, sample_name_getter ) ;
		m_have_written_header = true ;
		m_number_of_samples = number_of_samples ;
		m_streams.resize( m_number_of_samples ) ;
	}
	
	SNPDataSink::SinkPos VCFFormatSNPDataSink::get_stream_pos() const {
		if( m_compression_type == CompressionType( "no_compression" )) {
			return SinkPos( this, m_stream_ptr->tellp() ) ;
		}
		else {
			throw OperationUnsupportedError( "genfile::VCFFormatSNPDataSink::get_stream_pos()", "get stream position", "VCFFormatSNPDataSink( \"" + m_filename + "\" )" ) ;
		}
	}
}
