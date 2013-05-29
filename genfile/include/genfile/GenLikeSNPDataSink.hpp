
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENLIKESNPDATASINK_HPP
#define GENLIKESNPDATASINK_HPP

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSink.hpp"

namespace genfile {
	// This class represents a SNPDataSink which writes its data
	// to a flat file in a GEN-like format.
	// That is, SNP identifying data goes at the start of the line
	// followed by a sequence of numbers representing the data.
	class GenLikeSNPDataSink: public SNPDataSink
	{
	public:
		GenLikeSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type ) ;
		SinkPos get_stream_pos() const ;
		void omit_chromosome() ;
		std::string get_spec() const ;

	protected:
		
		operator bool() const ;
		std::ostream& stream() { return *m_stream_ptr ; }
		std::string const& filename() const { return m_filename ; }
		void write_variant( std::ostream& out, genfile::SNPIdentifyingData const& variant ) ;

	private:
		void setup( std::string const& filename, CompressionType compression_type ) ;

		std::string m_filename ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		Chromosome m_chromosome ;
		bool m_write_chromosome_column ;
		CompressionType m_compression_type ;
	} ;
}

#endif