#ifndef GENFILE_SORTINGBGENFILESNPDATASINK_HPP
#define GENFILE_SORTINGBGENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include <utility>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/GenomePosition.hpp"

namespace genfile {
	class SortingBGenFileSNPDataSink: public BasicBGenFileSNPDataSink
	{
	public:
		SortingBGenFileSNPDataSink(
			std::string const& filename,
			std::string const& free_data,
			bgen::uint32_t flags = bgen::e_CompressedSNPBlocks
		) ;

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

		~SortingBGenFileSNPDataSink() ;
	private:
		typedef std::string RSIDType ;
		typedef std::string DataType ;
		typedef std::pair< GenomePosition, RSIDType > SNPKey ;
		typedef std::pair< SNPKey, DataType > StoredRow ;
		std::vector< StoredRow > m_rows ;
	} ;
}

#endif
