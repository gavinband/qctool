
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_FORMAT_SNP_DATA_SINK_HPP
#define GENFILE_VCF_FORMAT_SNP_DATA_SINK_HPP

#include <string>
#include <set>
#include <iostream>
#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/SNPDataSink.hpp"

namespace genfile {
	struct VCFFormatSNPDataSink: public SNPDataSink {
	public:
		VCFFormatSNPDataSink( std::string const& filename ) ;

	public:		
		void set_output_fields( std::set< std::string > const& fields ) ;

	private:
		void write_header( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) const ;
		
		// Methods required by SNPDataSink
		operator bool() const { return *m_stream_ptr ; }
		
		void set_metadata_impl( Metadata const& metadata ) ;
		
		void write_variant_data_impl(
			VariantIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info
		) ;
		
		std::string get_spec() const ;

		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) ;
		SinkPos get_stream_pos() const ;

	private:
		std::string const m_filename ;
		boost::optional< std::set< std::string > > m_output_fields ;
		CompressionType const m_compression_type ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		bool m_have_written_header ;
		std::size_t m_number_of_samples ;
		boost::ptr_vector< std::ostringstream > m_streams ;
		Metadata m_metadata ;

	private:
		void write_info( Info const& info ) ;
	} ;
}

#endif
