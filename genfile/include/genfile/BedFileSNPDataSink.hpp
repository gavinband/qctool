
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_PEDIGREE_MAPPING_BEDFILE_SNP_DATA_SINK_HPP
#define GENFILE_PEDIGREE_MAPPING_BEDFILE_SNP_DATA_SINK_HPP

#include <vector>
#include <utility>
#include <map>
#include "genfile/SNPDataSink.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	
	/*
	* class PedgireeMappingBedFileSNPDataSink
	* Outputs data in PLINK Binary PED format (See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml)
	* Output is always in SNP-major mode.
	*/
	class BedFileSNPDataSink: public SNPDataSink
	{
	public:
		BedFileSNPDataSink(
			std::string const& output_filename,
			double call_threshhold = 0.9
		) ;

		// The destructor actually does the work of writing to the output file.
		~BedFileSNPDataSink() ;

		operator bool() const {
			return true ;
		}
		
		std::string get_spec() const ;

		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) ;
		void set_metadata_impl( Metadata const& ) ;

	private:
		std::vector< std::string > m_sample_ids ;
		std::string m_output_filename_stub ;
		double const m_call_threshhold ;
		
		std::auto_ptr< std::ostream > m_bim_file ;
		std::auto_ptr< std::ostream > m_bed_file ;
		
		std::vector< char > m_buffer ;
	private:
		void setup() ;
		void write_map_file( std::string const& output_filename ) const ;
		
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
			Info const& info
		) ;
	} ;
}

#endif
