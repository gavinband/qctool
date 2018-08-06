
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef LISTSNPDATASINK_HPP
#define LISTSNPDATASINK_HPP

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/Error.hpp"

namespace genfile {

	// This class encapsulates the basic method of writing a BGen file.
	class ListSNPDataSink: public SNPDataSink
	{
	public:
		// This class is intended to be used via a derived class.
		ListSNPDataSink(
			std::string const& filename,
			CompressionType compression_type
		) ;

		ListSNPDataSink(
			std::auto_ptr< std::ostream > stream_ptr,
			std::string const& filename,
			CompressionType compression_type
		) ;

		SinkPos get_stream_pos() const ;
		std::string get_spec() const ;

	public:
		// Methods required by SNPDataSink
		operator bool() const ;
		
		void write_variant_data_impl(
			VariantIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info
		) ;
		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) ;

	protected:
		// Other methods.
		std::auto_ptr< std::ostream >& stream_ptr() ;
		std::string const& filename() const ;
		
	private:

		void setup() ;

		std::string m_filename ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
	} ;
}

#endif

