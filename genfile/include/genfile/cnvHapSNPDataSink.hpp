
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CNVHAPSNPDATASINK_HPP
#define CNVHAPSNPDATASINK_HPP

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSink.hpp"

namespace genfile {
	// A SNPDataSink which writes its data in a format suitable for input
	// to cnvHap's normalisation pipeline, c.f.
	// http://www.nature.com/nmeth/journal/v7/n7/abs/nmeth.1466.html
	// http://www.imperial.ac.uk/people/l.coin
	//
	class cnvHapSNPDataSink: public SNPDataSink
	{
	public:
		cnvHapSNPDataSink( std::string const& filename ) ;
		SinkPos get_stream_pos() const ;
		std::string get_spec() const ;
		
	private:
		
		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) ;
		void set_metadata_impl( Metadata const& ) {} ;
		void write_variant_data_impl(
			SNPIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info
		) ;
		void finalise_impl() {}
		
		operator bool() const ;
		std::ostream& stream() { return *m_stream_ptr ; }
		std::string const& filename() const { return m_filename ; }

	private:
		void setup( std::string const& filename ) ;

		std::string m_filename ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		std::size_t m_number_of_samples ;
		Eigen::MatrixXd m_intensities ;
		Eigen::MatrixXd m_nonmissingness ;
	} ;
}

#endif
