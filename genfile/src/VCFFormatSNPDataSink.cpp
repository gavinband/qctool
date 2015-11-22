
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <boost/format.hpp>
#include <boost/bind.hpp>
#include "genfile/SNPDataSink.hpp"
#include "genfile/VCFFormatSNPDataSink.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"

#define DEBUG_VCFFORMATSNPDATASINK 1

namespace genfile {
	namespace {
		void get_first( std::vector< std::string >* target, std::string const& first, std::string const& ) {
			target->push_back( first ) ;
		}
	}
	
	VCFFormatSNPDataSink::VCFFormatSNPDataSink( std::string const& filename ):
		m_filename( filename ),
		m_compression_type( get_compression_type_indicated_by_filename( filename ) ),
		m_stream_ptr( open_text_file_for_output( filename, m_compression_type )),
		m_have_written_header( false ),
		m_number_of_samples( 0 )
	{
		(*m_stream_ptr) << std::resetiosflags( std::ios::floatfield ) << std::setprecision( 6 ) ;
#if 0
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "GT" ;
		format[ "Number" ] = "1" ;
		format[ "Description" ] = "Genotype call" ;
		m_metadata.insert( Metadata::value_type( "FORMAT", format )) ;
		format[ "ID" ] = "GP" ;
		format[ "Number" ] = "G" ;
		format[ "Description" ] = "Genotype call probabilities" ;
		m_metadata.insert( Metadata::value_type( "FORMAT", format )) ;
#endif
	}
	
	std::string VCFFormatSNPDataSink::get_spec() const {
		return m_filename ;
	}
	
	namespace {
		struct DataWriter: public VariantDataReader::PerSampleSetter {
			DataWriter( boost::ptr_vector< std::ostringstream >& streams, bool field_is_genotype ):
				m_streams( streams ),
				m_field_is_genotype( field_is_genotype ),
				m_sample_i( 0 ),
				m_sep( m_field_is_genotype ? '/' : ',' )
			{
			}
			
			~DataWriter() throw() {}

			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				assert( m_streams.size() == nSamples ) ;
				assert( nAlleles == 2 ) ;
				m_sample_i = 0 ;
			}

			bool set_sample( std::size_t i ) {
				m_sample_i = i ;
				return true ;
			}

			void set_number_of_entries( std::size_t n, OrderType const order_type, ValueType const value_type ) {
				m_number_of_entries = n ;
				m_order_type = order_type ;
				m_entry_i = 0 ;
				if( m_field_is_genotype ) {
					m_sep = ( order_type == ePerOrderedHaplotype ) ? '|' : '/' ;
				} else {
					m_sep = ',' ;
				}
			}

			void set_value( MissingValue const value ) {
				if( m_entry_i++ > 0 ) {
					m_streams[ m_sample_i ] << m_sep ;
				}
				m_streams[ m_sample_i ] << "." ;
			}

			void set_value( std::string& value ) {
				if( m_entry_i++ > 0 ) {
					m_streams[ m_sample_i ] << m_sep ;
				}
				m_streams[ m_sample_i ] << value ;
			}

			void set_value( Integer const value ) {
				if( m_entry_i++ > 0 ) {
					m_streams[ m_sample_i ] << m_sep ;
				}
				m_streams[ m_sample_i ] << value ;
			}

			void set_value( double const value ) {
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

	void VCFFormatSNPDataSink::set_metadata_impl( Metadata const& metadata ) {
		m_metadata = metadata ;
	}

	void VCFFormatSNPDataSink::set_output_fields( std::set< std::string > const& fields ) {
		m_output_fields = fields ;
	}

	namespace impl {
		void write_metadata(
			boost::optional< std::set< std::string > > output_fields,
			SNPDataSink::Metadata const& metadata,
			std::ostream& stream
		) {
			
		}
	}

	void VCFFormatSNPDataSink::write_header( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) const {
		(*m_stream_ptr) << "##fileformat=VCFv4.1\n" ;
		
		std::pair< Metadata::const_iterator, Metadata::const_iterator > const
			formatDefinitions = std::make_pair( m_metadata.begin(), m_metadata.end() ) ;

		boost::format formatFormat( "##FORMAT=<ID=%s,Type=%s,Number=%s,Description=\"%s\">\n" ) ;

#if DEBUG_VCFFORMATSNPDATASINK
		Metadata::const_iterator format_i = m_metadata.begin() ;
		std::cerr << "VCFFormatSNPDataSink::write_header(): FORMAT entries are:\n" ;
		for( ; format_i != m_metadata.end(); ++format_i ) {
			if( format_i->first == "FORMAT" ) {
				std::cerr << ( formatFormat % format_i->second.at( "ID" ) % format_i->second.at( "Type" ) % format_i->second.at( "Number" ) % format_i->second.at( "Description" ) ) ;
			}
		}
		std::cerr << "\n" ;
#endif

		if( !m_output_fields ) {
			Metadata::const_iterator format_i = formatDefinitions.first ;
			for( ; format_i != formatDefinitions.second; ++format_i ) {
				if( format_i->first == "FORMAT" ) {
					(*m_stream_ptr) << ( formatFormat % format_i->second.at( "ID" ) % format_i->second.at( "Type" ) % format_i->second.at( "Number" ) % format_i->second.at( "Description" ) ) ;
				}
			}
		} else {
			std::set< std::string >::const_iterator i = m_output_fields->begin() ;
			std::set< std::string >::const_iterator end_i = m_output_fields->end() ;
			for( ; i != end_i; ++i ) {
				Metadata::const_iterator format_i = formatDefinitions.first ;
				bool found = false ;
				for( ; format_i != formatDefinitions.second; ++format_i ) {
					if( format_i->first == "FORMAT" && format_i->second.at( "ID" ) == *i ) {
						(*m_stream_ptr) << ( formatFormat % format_i->second.at( "ID" ) % format_i->second.at( "Type" ) % format_i->second.at( "Number" ) % format_i->second.at( "Description" ) ) ;
						found = true ;
					}
				}
				if( !found ) {
					std::map< std::string, std::string > unknownFormat ;
					unknownFormat[ "ID" ] = *i ;
					unknownFormat[ "Number" ] = "." ;
					unknownFormat[ "Description" ] = "Unknown field" ;
					(*m_stream_ptr) << ( formatFormat % *i % "String" % "." % "\"Unknown field\"" ) ;
				}
			}
		}

		(*m_stream_ptr) << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			(*m_stream_ptr) << "\t" << sample_name_getter( i ) ;
		}
		(*m_stream_ptr) << "\n" ;
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
			<< id_data.get_rsid() ;
		if( id_data.get_SNPID() != id_data.get_rsid() ) {
			(*m_stream_ptr ) << "," << id_data.get_SNPID() ;
		}
		(*m_stream_ptr) << tab
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
			if(
				( !(m_output_fields) || ( m_output_fields.get().find( field ) != m_output_fields.get().end() ) )
				&& field != ":genotypes:"
				&& field != ":intensities:"
			) {
				if( data_reader.supports( field )) {
					(*m_stream_ptr) << ( have_format ? ":" : "" ) << field ;
					DataWriter writer( m_streams, ( field == "GT" ) ) ; // kludge, fix this.
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
