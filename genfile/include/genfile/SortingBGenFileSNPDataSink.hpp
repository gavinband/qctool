
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SORTINGBGENFILESNPDATASINK_HPP
#define GENFILE_SORTINGBGENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include <utility>
#include <map>
#include <stdint.h>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/GenomePosition.hpp"

namespace genfile {
	// This class constructs a Bgen file which writes its data to a temporary file,
	// and then copies that file in the right order into the destination.
	class SortingBGenFileSNPDataSink: public SNPDataSink
	{
	public:
		SortingBGenFileSNPDataSink(
			std::string const& filename, 
			SNPDataSink::UniquePtr m_sink
		) ;

		void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability,
			Info const&
		) ;

		~SortingBGenFileSNPDataSink() ;
		
		std::string get_spec() const ;
		
		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) ;
		
	public:
		// return the number of samples represented in SNPs in the file.
		// The value returned is undefined until after the first snp has been written.
		uint32_t number_of_samples() const { return m_sink->number_of_samples() ; }
		// return the number of SNPs that have been written to the file so far.
		std::size_t number_of_snps_written() const { return m_sink->number_of_snps_written() ; }

	public:
		// The following functions must be implemented by derived classes.
		operator bool() const { return m_sink.get() && (*m_sink) ; }
		
	private:
		std::string m_filename ;
		SNPDataSink::UniquePtr m_sink ;
		typedef std::multimap< VariantIdentifyingData, std::pair< std::ostream::streampos, std::ostream::streampos > > OffsetMap ;
		OffsetMap m_file_offsets ;
		std::ostream::streampos m_offset_of_first_snp ;
	} ;
}

#endif
