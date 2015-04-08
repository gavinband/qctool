
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
	class PedigreeMappingPedFileSNPDataSink: public SNPDataSink
	{
	public:
		PedigreeMappingPedFileSNPDataSink(
			CohortIndividualSource const& samples,
			Pedigree const& pedigree,
			std::string const& output_filename,
			double call_threshhold = 0.9
		) ;

		// The destructor actually does the work of writing to the output file.
		~PedigreeMappingPedFileSNPDataSink() ;

		operator bool() const {
			return true ;
		}

		std::string get_spec() const ;
	private:
		CohortIndividualSource const& m_samples ;
		Pedigree const& m_pedigree ;
		std::vector< std::string > const m_phenotypes ;
		std::map< std::string, std::size_t > m_pedigree_to_sample_mapping ;
		std::string m_output_filename_stub ;
		std::vector< SNPIdentifyingData > m_written_snps ;
		std::vector< std::vector< std::pair< char, char > > > m_written_alleles ;
		double const m_call_threshhold ;
	private:
		static std::map< std::string, std::size_t > get_pedigree_to_sample_mapping( Pedigree const& pedigree, CohortIndividualSource const& samples ) ;

		// Write PED file suitable for PLINK or QTDT
		void write_ped_file( std::string const& output_filename ) const ;
		// Write DAT file suitable for QTDT
		void write_dat_file( std::string const& output_filename ) const ;
		// Write MAP file suitable for plink --map3
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
