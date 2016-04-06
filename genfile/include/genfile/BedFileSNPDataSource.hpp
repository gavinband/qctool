
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_BEDFILESNPDATAPROVIDER_HPP
#define GENFILE_BEDFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "IdentifyingDataCachingSNPDataSource.hpp"

namespace genfile {
	namespace impl {
		struct BedFileSNPDataReader ;
	}

	// This class represents a SNPDataSource which reads its data
	// from a BGEN file.
	class BedFileSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
		friend struct impl::BedFileSNPDataReader ;
	public:
		BedFileSNPDataSource( std::string const& bedFilename, std::string const& bimFilename, std::string const& famFilename ) ;

		Metadata get_metadata() const ;
		unsigned int number_of_samples() const { return m_number_of_samples ; }
		OptionalSnpCount total_number_of_snps() const { return m_total_number_of_snps ; }
		operator bool() const { return !m_exhausted && *m_bed_stream_ptr ; }

		std::istream& stream() { return *m_bed_stream_ptr ; }
		std::istream const& stream() const { return *m_bed_stream_ptr ; }
		std::string get_source_spec() const { return m_bed_filename ; }

	private:

		void reset_to_start_impl() ;

		void read_snp_identifying_data_impl( VariantIdentifyingData* result ) ;
		VariantDataReader::UniquePtr read_variant_data_impl() ;
		void ignore_snp_probability_data_impl() ;

	private:
		std::string const m_bed_filename ;
		std::string const m_bim_filename ;
		bool m_exhausted ;
		unsigned int m_number_of_samples ;
		OptionalSnpCount m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_bed_stream_ptr ;
		std::auto_ptr< std::istream > m_bim_stream_ptr ;
		std::vector< std::pair< int64_t, int64_t > > m_genotype_table ;
		void setup( std::string const& bedFilename, std::string const& bimFilename, std::string const& famFilename ) ;
	} ;
}

#endif
