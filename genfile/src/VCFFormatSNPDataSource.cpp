#include <string>
#include <map>
#include <memory>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/VCFFormatSNPDataSource.hpp"

using boost::tuples::tie ;

namespace genfile {
	VCFFormatSNPDataSource::VCFFormatSNPDataSource( std::auto_ptr< std::istream > stream_ptr ):
		m_spec( "(unnamed stream)" ),
		m_stream_ptr( stream_ptr )
	{
		assert( *m_stream_ptr ) ;
	}

	VCFFormatSNPDataSource::VCFFormatSNPDataSource( std::string const& filename ):
		m_spec( "file://" + filename ),
		m_stream_ptr( open_text_file_for_input( filename, get_compression_type_indicated_by_filename( filename )))
	{}

	VCFFormatSNPDataSource::operator bool() const {
		return *m_stream_ptr ;
	}

	unsigned int VCFFormatSNPDataSource::number_of_samples() const {
		assert(0) ;
	}

	unsigned int VCFFormatSNPDataSource::total_number_of_snps() const {
		assert(0) ;
	}

	std::string VCFFormatSNPDataSource::get_source_spec() const {
		return m_spec ;
	}

	std::string VCFFormatSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + m_spec ;
	}

	void VCFFormatSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		assert(0) ;
	}

	void VCFFormatSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		assert(0) ;
	}

	void VCFFormatSNPDataSource::ignore_snp_probability_data_impl() {
		assert(0) ;
	}

	void VCFFormatSNPDataSource::reset_to_start_impl() {
		assert(0) ;
	}
}
