
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BIMBAMDOSAGESNPDATASINK_HPP
#define BIMBAMDOSAGESNPDATASINK_HPP

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/GenLikeSNPDataSink.hpp"

namespace genfile {
	
	// This class represents a SNPDataSink which writes its data
	// to a plain GEN file.
	class BimBamDosageSNPDataSink: public GenLikeSNPDataSink
	{
	public:
		BimBamDosageSNPDataSink( std::string const& filename ) ;
		BimBamDosageSNPDataSink( std::string const& filename, CompressionType compression_type ) ;

		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) ;
		void write_variant( std::ostream& out, genfile::VariantIdentifyingData const& variant ) ; // overrides
	private:

		Eigen::VectorXd m_data ;
		
		void write_variant_data_impl(
			VariantIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info = Info()
		) ;
	} ;
}

#endif
