#ifndef GENFILE_PEDFILE_SNP_DATA_SINK_HPP
#define GENFILE_PEDFILE_SNP_DATA_SINK_HPP

#include <vector>
#include <utility>
#include <map>
#include "genfile/SNPDataSink.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/Pedigree.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {
	class PedFileSNPDataSink: public SNPDataSink
	{
	public:
		PedFileSNPDataSink(
			CohortIndividualSource const& samples,
			Pedigree const& pedigree,
			std::string const& phenotype,
			std::string const& output_filename,
			double call_threshhold = 0.9
		) ;

		// The destructor actually does the work of writing to the output file.
		~PedFileSNPDataSink() ;

		operator bool() const {
			return true ;
		}
		
	private:
		CohortIndividualSource const& m_samples ;
		Pedigree const& m_pedigree ;
		std::string const m_phenotype ;
		std::map< std::string, std::size_t > m_pedigree_to_sample_mapping ;
		std::string const m_output_filename ;
		std::vector< SNPIdentifyingData > m_written_snps ;
		std::vector< std::vector< std::pair< char, char > > > m_written_alleles ;
		double const m_call_threshhold ;
	private:
		
		static std::map< std::string, std::size_t > get_pedigree_to_sample_mapping( Pedigree const& pedigree, CohortIndividualSource const& samples ) ;
		void write_ped_file( std::string const& output_filename ) const ;
		void write_map_file( std::string const& output_filename ) const ;
		
		void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) ;
	} ;
}

#endif
