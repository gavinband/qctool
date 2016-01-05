
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef PENNCNVSNPDATASINK_HPP
#define PENNCNVSNPDATASINK_HPP

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
	class PennCNVSNPDataSink: public SNPDataSink
	{
	public:
		PennCNVSNPDataSink( std::string const& filename, double const threshhold = 0.9 ) ;
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
		double const m_threshhold ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		Eigen::MatrixXd m_intensities ;
		Eigen::MatrixXd m_nonmissingness ;
		std::vector< std::size_t > m_genotypes ;
		std::vector< std::string > m_genotype_alleles ;
	} ;
}

#endif
