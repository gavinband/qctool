
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/SNPDataSource.hpp"
#include "genfile/Error.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/GenomePosition.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "statfile/SNPDataSourceAdapter.hpp"

namespace statfile {
	SNPDataSourceAdapter::SNPDataSourceAdapter( BuiltInTypeStatSource::UniquePtr source ):
		m_source( source )
	{
		setup() ;
	}

	void SNPDataSourceAdapter::setup() {
		using genfile::string_utils::to_lower ;
		assert( m_source.get() ) ;
		if( m_source->number_of_columns() < 6 ) {
			throw genfile::MalformedInputError( get_source_spec(), 0 ) ;
		}
		// Could check all the columns.  But let's allow some flexibility.
		if( to_lower( m_source->name_of_column(0) ) != "snpid" ) {
			throw genfile::MalformedInputError( get_source_spec(), 0, 0 ) ;
		}
		if( to_lower( m_source->name_of_column(1) ) != "rsid" ) {
			throw genfile::MalformedInputError( get_source_spec(), 0, 1 ) ;
		}
		if( to_lower( m_source->name_of_column(2) ) != "chromosome" && to_lower( m_source->name_of_column(2) ) != "chr" ) {
			throw genfile::MalformedInputError( get_source_spec(), 0, 2 ) ;
		}
		if( to_lower( m_source->name_of_column(3) ) != "position" && to_lower( m_source->name_of_column(3) ) != "pos" ) {
			throw genfile::MalformedInputError( get_source_spec(), 0, 3 ) ;
		}
	}

	std::string SNPDataSourceAdapter::get_source_spec() const {
		return "SNPDataSourceAdapter( " + m_source->get_source_spec() + " )" ;
	}

	std::string SNPDataSourceAdapter::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "SNPDataSourceAdapter( " + m_source->get_source_spec() + " )" ;
	}

	void SNPDataSourceAdapter::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		if( m_source->current_column() == 0 ) {
			std::string SNPID, rsid, chr_string, alleleA, alleleB ;
			uint32_t position ;
			(*m_source) >> SNPID >> rsid >> chr_string >> position >> alleleA >> alleleB ;
			if( *this ) {
				m_snp = genfile::SNPIdentifyingData( SNPID, rsid, genfile::GenomePosition( chr_string, position ), alleleA, alleleB ) ;
			}
		}
		if( *this ) {
			set_number_of_samples( 0 ) ;
			set_SNPID( m_snp.get_SNPID() ) ;
			set_RSID( m_snp.get_rsid() ) ;
			set_chromosome( m_snp.get_position().chromosome() ) ;
			set_SNP_position( m_snp.get_position().position() ) ;
			set_allele1( m_snp.get_first_allele() ) ;
			set_allele2( m_snp.get_second_allele() ) ;
		}
	}

	void SNPDataSourceAdapter::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& ,
		std::string const& genotype_field
	) {
		assert( m_source->current_column() == 6 ) ;
		(*m_source) >> ignore_all() ;
	}

	namespace {
		struct DataReader: public genfile::VariantDataReader {
			VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				setter.set_number_of_samples( 0 ) ;
				return *this ;
			}
			VariantDataReader& get( std::string const& spec, PerVariantSetter& setter ) {
				return *this ;
			}
			bool supports( std::string const& spec ) const {
				return false ;
			}
			void get_supported_specs( SpecSetter ) const {
				return ;
			}
			std::size_t get_number_of_samples() const {
				return 0 ;
			}
		} ;
	}

	genfile::VariantDataReader::UniquePtr SNPDataSourceAdapter::read_variant_data_impl() {
		return genfile::VariantDataReader::UniquePtr( new DataReader() ) ;
	}

	void SNPDataSourceAdapter::ignore_snp_probability_data_impl() {
		assert( m_source->current_column() == 6 ) ;
		(*m_source) >> ignore_all() ;
	}

	void SNPDataSourceAdapter::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}

