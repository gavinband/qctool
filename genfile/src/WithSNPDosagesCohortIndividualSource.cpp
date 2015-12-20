
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <string>
#include <map>
#include <iomanip>


#include "../config.hpp"
#if HAVE_EIGEN
#include <Eigen/Core>

#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CompositeCohortIndividualSource.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/get_set.hpp"
#include "genfile/get_set_eigen.hpp"
#include "genfile/WithSNPDosagesCohortIndividualSource.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"
#include "genfile/SNPIDMatchesTest.hpp"


namespace genfile {
	VariantIdentifyingDataTest::UniquePtr WithSNPDosagesCohortIndividualSource::create_snp_matcher( std::string test_spec, std::string const& splitter ) {
		std::vector< std::string >  bits = genfile::string_utils::split_and_strip( test_spec, splitter, " " ) ;
		if( bits.size() == 0 || bits.size() > 2 ) {
			throw genfile::BadArgumentError( "SNPMatcher::create()", "test_spec = \"" + test_spec + "\"" ) ;
		}
		// If no ~, default to rsid matcher.
		if( bits.size() == 1 ) {
			bits.insert( bits.begin(), "rsid" ) ;
		}

		assert( bits.size() == 2 ) ;

		if( bits[0] == "rsid" ) {
			return VariantIdentifyingDataTest::UniquePtr( new SNPIDMatchesTest( "rsid~" + bits[1] )) ;
		}
		else if( bits[0] == "snpid" ) {
			return VariantIdentifyingDataTest::UniquePtr( new SNPIDMatchesTest( "snpid~" + bits[1] )) ;
		}
		else if( bits[0] == "pos" || bits[0] == "position" ) {
			return VariantIdentifyingDataTest::UniquePtr( new PositionMatchesTest( bits[1] )) ;
		}
		else {
			throw genfile::BadArgumentError( "genfile::create_matcher()", "test_spec = \"" + test_spec + "\"" ) ;
		}
	}
	
	PositionMatchesTest::PositionMatchesTest( GenomePosition const& position ):
		m_position( position )
	{}

	bool PositionMatchesTest::operator()( VariantIdentifyingData const& data ) const {
		return data.get_position() == m_position ;
	}
		
	std::string PositionMatchesTest::display() const {
		return "position=" + genfile::string_utils::to_string( m_position ) ;
	}

	WithSNPDosagesCohortIndividualSource::UniquePtr WithSNPDosagesCohortIndividualSource::create(
		genfile::CohortIndividualSource::ConstUniquePtr sample_source,
		genfile::SNPDataSource& snp_data_source,
		SNPDosageSpec const& snp_matchers
	) {
		return UniquePtr(
			new WithSNPDosagesCohortIndividualSource(
				sample_source,
				snp_data_source,
				snp_matchers
			)
		) ;
	}
		
	WithSNPDosagesCohortIndividualSource::WithSNPDosagesCohortIndividualSource(
		genfile::CohortIndividualSource::ConstUniquePtr sample_source,
		genfile::SNPDataSource& snp_data_source,
		SNPDosageSpec const& snp_matchers
	):
		m_source( sample_source )
	{
		assert( m_source->get_number_of_individuals() == snp_data_source.number_of_samples() ) ;
		setup( snp_data_source, snp_matchers ) ;
	}
		
	std::size_t WithSNPDosagesCohortIndividualSource::get_number_of_individuals() const {
		return m_source->get_number_of_individuals() ;
	}
	
	genfile::CohortIndividualSource::ColumnSpec WithSNPDosagesCohortIndividualSource::get_column_spec() const {
		return m_source->get_column_spec() + ColumnSpec( m_column_names, m_column_types ) ;
	}
	
	genfile::CohortIndividualSource::Entry WithSNPDosagesCohortIndividualSource::get_entry( std::size_t sample_i, std::string const& column_name ) const {
		assert( sample_i < get_number_of_individuals() ) ;
		if( m_source->check_for_column( column_name )) {
			return m_source->get_entry( sample_i, column_name ) ;
		}
		else {
			std::vector< std::string >::const_iterator where = std::find( m_column_names.begin(), m_column_names.end(), column_name ) ;
			assert( where != m_column_names.end() ) ;
			return m_column_data[ std::size_t( where - m_column_names.begin() ) ][ sample_i ] ;
		}
	}
	
	genfile::CohortIndividualSource const& WithSNPDosagesCohortIndividualSource::get_parent_source() const {
		return *m_source ;
	}

	genfile::CohortIndividualSource const& WithSNPDosagesCohortIndividualSource::get_base_source() const {
		return m_source->get_base_source() ;
	}
	
	std::string WithSNPDosagesCohortIndividualSource::get_source_spec() const {
		return "with-snp-dosages:" + m_source->get_source_spec() ;
	}
	
	void WithSNPDosagesCohortIndividualSource::setup(
		genfile::SNPDataSource& snp_data_source,
		SNPDosageSpec const& snp_matchers
	) {
		std::map< std::string, std::size_t > snp_counts ;
		VariantIdentifyingData snp ;
		while(
			snp_data_source.get_snp_identifying_data( &snp )
		) {
			SNPDosageSpec::const_iterator where = snp_matchers.begin() ;
			for( ; where != snp_matchers.end(); ++where ) {
				if( where->first->operator()( snp )) {
					break ;
				}
			}
			if( where != snp_matchers.end() ) {
				Eigen::MatrixXd probabilities( snp_data_source.number_of_samples(), 3 ) ;
				if( !snp_data_source.read_snp_probability_data(
						genfile::set_genotypes( probabilities )
				)) {
					throw genfile::MalformedInputError(
						snp_data_source.get_source_spec(),
						snp_data_source.number_of_snps_read()
					) ;
				}
				
				Eigen::VectorXd row_sums = probabilities.rowwise().sum() ;
				for( int i = 0; i < row_sums.size(); ++i ) {
					if( row_sums(i) == 0.0 ) {
						probabilities.row( i ) = Eigen::Vector3d::Constant( std::numeric_limits< double >::quiet_NaN() ) ;
					}
				}
				
				std::string const key = where->first->display() ;
				
				if( ( ++( snp_counts[ key ] )) > 1 ) {
					throw genfile::DuplicateSNPError( snp_data_source.get_source_spec(), key ) ;
				}

				if( where->second.find( "add" ) != where->second.end() ) {
					add_column(
						key + ":additive_dosage",
						( probabilities.col(1) + 2.0 * probabilities.col(2) )
					) ;
				}

				if( where->second.find( "dom" ) != where->second.end() ) {
					add_column(
						key + ":dominant_dosage",
						( probabilities.col(1) + probabilities.col(2) )
					) ;
				}

				if( where->second.find( "het" ) != where->second.end() ) {
					add_column(
						key + ":heterozygote_dosage",
						( probabilities.col(1) )
					) ;
				}

				if( where->second.find( "rec" ) != where->second.end() ) {
					add_column(
						key + ":recessive_dosage",
						( probabilities.col(2) )
					) ;
				}
						}
			else {
				snp_data_source.ignore_snp_probability_data() ;
			}
		}

		// Check SNP counts are all greater than zero.
		for( SNPDosageSpec::const_iterator where = snp_matchers.begin() ; where != snp_matchers.end(); ++where ) {
			if( snp_counts[ where->first->display() ] == 0 ) {
				throw genfile::KeyNotFoundError( snp_data_source.get_source_spec(), where->first->display() ) ;
			}
		}
	}
	
	void WithSNPDosagesCohortIndividualSource::add_column( std::string const& column_name, Eigen::MatrixXd const& column ) {
		if( check_for_column( column_name ) ) {
			throw genfile::ColumnAlreadyExistsError( get_source_spec(), column_name ) ;
		}
		assert( std::size_t( column.rows() ) == get_number_of_individuals() ) ;
		m_column_names.push_back( column_name ) ;
		m_column_types.push_back( e_CONTINUOUS_COVARIATE ) ;
		m_column_data.push_back( std::vector< Entry >( column.data(), column.data() + column.rows() )) ;
		for( int i = 0;i < column.rows(); ++i ) {
			double value = m_column_data.back()[i].as< double >() ;
			if( value != value ) {
				m_column_data.back()[i] = MissingValue() ;
			}
		}
	}
}

#endif
