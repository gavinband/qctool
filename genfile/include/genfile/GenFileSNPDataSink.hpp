
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILESNPDATASINK_HPP
#define GENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/GenLikeSNPDataSink.hpp"

namespace genfile {
	
	// This class represents a SNPDataSink which writes its data
	// to a plain GEN file.
	class GenFileSNPDataSink: public GenLikeSNPDataSink
	{
	public:
		GenFileSNPDataSink( std::string const& filename, Chromosome chromosome ) ;
		GenFileSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type ) ;

		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) ;

	private:
		Eigen::VectorXd m_data ;
		
		void write_variant_data_impl(
			SNPIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info = Info()
		) ;
	} ;
}

#endif
