
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <boost/bind.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/Error.hpp"
#include "genfile/ListSNPDataSink.hpp"
#include "genfile/ToGP.hpp"

namespace genfile {
	// This class is intended to be used via a derived class.
	ListSNPDataSink::ListSNPDataSink(
		std::string const& filename,
		CompressionType compression_type
	)
	: 	m_filename( filename ),
		m_stream_ptr( open_text_file_for_output( m_filename, compression_type ) )
	{
		setup() ;
	}

	ListSNPDataSink::ListSNPDataSink(
		std::auto_ptr< std::ostream > stream_ptr,
		std::string const& filename,
		CompressionType compression_type
	)
	: 	m_filename( filename ),
		m_stream_ptr( stream_ptr )
	{
		setup() ;
	}

	void ListSNPDataSink::setup() {
		(*m_stream_ptr) << "SNPID rsid chromosome position alleleA alleleB\n" ;
	}

	std::string ListSNPDataSink::get_spec() const {
		return "list(\"" + m_filename + "\")" ;
	}

	void ListSNPDataSink::set_sample_names_impl( std::size_t, SampleNameGetter ) {
	}

	SNPDataSink::SinkPos ListSNPDataSink::get_stream_pos() const {
		return SinkPos( this, m_stream_ptr->tellp() ) ;
	}
	
	ListSNPDataSink::operator bool() const { return m_stream_ptr->good() ; }

	void ListSNPDataSink::write_variant_data_impl(
		VariantIdentifyingData const& id_data,
		VariantDataReader&,
		Info const&
	) {
		std::string const& SNPID = ( id_data.number_of_identifiers() > 1 ? id_data.get_identifiers_as_string(",", 1) : id_data.get_identifiers_as_string(",", 0,1) ) ;
		std::string chromosome = "NA" ;
		if( !id_data.get_position().chromosome().is_missing() ) {
			chromosome = id_data.get_position().chromosome() ;
		}
		(*m_stream_ptr) << SNPID
			<< " " << id_data.get_primary_id()
			<< " " << chromosome
			<< " " << id_data.get_position().position()
			<< " " << id_data.get_allele(0)
			<< " " << id_data.get_alleles_as_string(",",1)
			<< "\n" ;
	}

}

