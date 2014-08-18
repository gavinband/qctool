
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
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
		
		void get_first( std::vector< std::string >* target, std::string const& first, std::string const& ) {
			target->push_back( first ) ;
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
			DataWriter( boost::ptr_vector< std::ostringstream >& streams, bool field_is_genotype ):
				m_streams( streams ),
				m_field_is_genotype( field_is_genotype ),
				m_sep( m_field_is_genotype ? '/' : ',' ),
				m_sample_i( 0 )
			{
			}
			
			~DataWriter() throw() {}

			void set_number_of_samples( std::size_t n ) {
				assert( m_streams.size() == n ) ;
				m_sample_i = 0 ;
			}

			void set_sample( std::size_t i ) {
				m_sample_i = i ;
			}

			void set_number_of_entries( std::size_t n ) {
				m_number_of_entries = n ;
				m_entry_i = 0 ;
			}

			void set_order_type( OrderType const type ) {
				m_order_type = type ;
				if( m_field_is_genotype ) {
					m_sep = ( type == eOrderedList ) ? '|' : '/' ;
				} else {
					m_sep = ',' ;
				}
			}

			void operator()( MissingValue const value ) {
				if( m_entry_i++ > 0 ) {
					m_streams[ m_sample_i ] << m_sep ;
				}
				m_streams[ m_sample_i ] << "." ;
			}

			void operator()( std::string& value ) {
				if( m_entry_i++ > 0 ) {
					m_streams[ m_sample_i ] << m_sep ;
				}
				m_streams[ m_sample_i ] << value ;
			}

			void operator()( Integer const value ) {
				if( m_entry_i++ > 0 ) {
					m_streams[ m_sample_i ] << m_sep ;
				}
				m_streams[ m_sample_i ] << value ;
			}

			void operator()( double const value ) {
				if( m_entry_i++ > 0 ) {
					m_streams[ m_sample_i ] << m_sep ;
				}
				m_streams[ m_sample_i ] << value ;
			}

		private:
			boost::ptr_vector< std::ostringstream >& m_streams ;
			bool m_field_is_genotype ;
			std::size_t m_sample_i ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
			OrderType m_order_type ;
			char m_sep ;
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

		std::vector< std::string > fields ;
		data_reader.get_supported_specs( boost::bind( get_first, &fields, _1, _2 ) ) ;
		bool have_format = false ;
		(*m_stream_ptr) << tab ;
		for( std::size_t field_i = 0; field_i < fields.size(); ++field_i ) {
			std::string const field = fields[ field_i ] ;
			if( field != ":genotypes:" && field != ":intensities:" ) {
				if( data_reader.supports( field )) {
					(*m_stream_ptr) << ( have_format ? ":" : "" ) << field ;
					DataWriter writer( m_streams, ( field == "GT" ) ) ;
					if( have_format ) {
						for( std::size_t i = 0; i < m_streams.size(); ++i ) {
							m_streams[i] << ":" ;
						}
					}
					data_reader.get( fields[field_i], writer ) ;
					have_format = true ;
				}
			}
		}

		for( std::size_t i = 0; i < m_streams.size(); ++i ) {
			if( m_streams[i].str() == "" ) {
				(*m_stream_ptr) << tab << "." ;
			} else {
				(*m_stream_ptr) << tab << m_streams[i].str() ;
			}
			m_streams[i].str( "" ) ;
		}
		(*m_stream_ptr) << "\n" ;
	}

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
