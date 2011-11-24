#ifndef GENFILE_BEDFILE_SNP_DATA_SINK_HPP
#define GENFILE_BEDFILE_SNP_DATA_SINK_HPP

#include <vector>
#include <utility>
#include <map>
#include "genfile/SNPDataSink.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/Pedigree.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	
	/*
	* class BedFileSNPDataSink
	* Outputs data in PLINK Binary PED format (See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml)
	* Output is always in SNP-major mode.
	*/
	class BedFileSNPDataSink: public SNPDataSink
	{
	public:
		BedFileSNPDataSink(
			CohortIndividualSource const& samples,
			Pedigree const& pedigree,
			std::string const& output_filename,
			double call_threshhold = 0.9
		) ;

		// The destructor actually does the work of writing to the output file.
		~BedFileSNPDataSink() ;

		operator bool() const {
			return true ;
		}

	private:
		CohortIndividualSource const& m_samples ;
		Pedigree const& m_pedigree ;
		std::vector< std::string > const m_phenotypes ;
		std::map< std::size_t, std::size_t > m_pedigree_to_sample_mapping ;
		std::string m_output_filename_stub ;
		double const m_call_threshhold ;
		
		std::auto_ptr< std::ostream > m_bim_file ;
		std::auto_ptr< std::ostream > m_bed_file ;
	private:
		void setup() ;
		std::map< std::size_t, std::size_t > get_pedigree_to_sample_mapping( Pedigree const& pedigree, CohortIndividualSource const& samples ) ;
		void write_ped_file( std::string const& output_filename ) const ;
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
			GenotypeProbabilityGetter const& get_BB_probability
		) ;
	} ;
}

#endif
