
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <set>
#include <cassert>
#include <exception>
#include <iostream>
#include <fstream>
#include <boost/variant.hpp>

#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/MissingValue.hpp"
#include "genfile/FromFileCohortIndividualSource.hpp"
#include "genfile/TraditionalStrictCohortIndividualSource.hpp"
#include "genfile/CategoricalCohortIndividualSource.hpp"

namespace genfile {
	CohortIndividualSource::SingleColumnSpec::SingleColumnSpec( std::string const& name, CohortIndividualSource::ColumnType const& type ):
		Base( name, type )
	{}

	CohortIndividualSource::SingleColumnSpec::SingleColumnSpec( CohortIndividualSource::SingleColumnSpec const& other ):
		Base( other )
	{}
	
	CohortIndividualSource::SingleColumnSpec& CohortIndividualSource::SingleColumnSpec::operator=( CohortIndividualSource::SingleColumnSpec const& other ) {
		Base::operator=( other ) ;
		return *this ;
	}

	std::string const& CohortIndividualSource::SingleColumnSpec::name() const { return first ; }

	CohortIndividualSource::ColumnType const CohortIndividualSource::SingleColumnSpec::type() const { return second ; }

	bool CohortIndividualSource::SingleColumnSpec::is_discrete() const {
		return second == e_DISCRETE_COVARIATE || second == e_BINARY_PHENOTYPE ;
	}

	bool CohortIndividualSource::SingleColumnSpec::is_continuous() const {
		return second == e_CONTINUOUS_COVARIATE || second == e_CONTINUOUS_PHENOTYPE ;
	}
	
	bool CohortIndividualSource::SingleColumnSpec::is_phenotype() const {
		return second == e_BINARY_PHENOTYPE || second == e_CONTINUOUS_PHENOTYPE ;
	}
	
	bool CohortIndividualSource::SingleColumnSpec::is_covariate() const {
		return second == e_DISCRETE_COVARIATE || second == e_CONTINUOUS_COVARIATE ;
	}
	
	
	bool CohortIndividualSource::SingleColumnSpec::operator==( CohortIndividualSource::SingleColumnSpec const& right ) const {
		return Base( *this ) == Base( right ) ;
	}
	
	bool CohortIndividualSource::SingleColumnSpec::operator!=( CohortIndividualSource::SingleColumnSpec const& right ) const {
		return Base( *this ) != Base( right ) ;
	}
	
	CohortIndividualSource::ColumnSpec::ColumnSpec() {}	

	CohortIndividualSource::ColumnSpec::ColumnSpec( std::vector< std::string > const& column_names, std::vector< ColumnType > const& column_types ):
		m_column_names( column_names ),
		m_column_types( column_types )
	{
		assert( m_column_names.size() == m_column_types.size() ) ;
	}

	CohortIndividualSource::ColumnSpec::ColumnSpec( ColumnSpec const& other ):
		m_column_names( other.m_column_names ),
		m_column_types( other.m_column_types )
	{}
	

	CohortIndividualSource::ColumnSpec& CohortIndividualSource::ColumnSpec::operator=( ColumnSpec const& other ) {
		m_column_names = other.m_column_names ;
		m_column_types = other.m_column_types ;
		return *this ;
	}

	CohortIndividualSource::SingleColumnSpec CohortIndividualSource::ColumnSpec::get_spec( std::size_t i ) const {
		assert( i < m_column_names.size() ) ;
		return SingleColumnSpec( m_column_names[i], m_column_types[i] ) ;
	}


	CohortIndividualSource::SingleColumnSpec CohortIndividualSource::ColumnSpec::operator[]( std::size_t i ) const {
		return get_spec( i ) ;
	}

	CohortIndividualSource::SingleColumnSpec CohortIndividualSource::ColumnSpec::operator[]( std::string const& column_name ) const {
		return get_spec( find_column( column_name ) ) ;
	}

	std::vector< std::string > CohortIndividualSource::ColumnSpec::get_names() const {
		return m_column_names ;
	}
	
	std::vector< CohortIndividualSource::ColumnType > CohortIndividualSource::ColumnSpec::get_types() const {
		return m_column_types ;
	}

	std::size_t CohortIndividualSource::ColumnSpec::size() const {
		return m_column_names.size() ;
	}

	bool CohortIndividualSource::ColumnSpec::operator==( ColumnSpec const& other ) {
		return m_column_names == other.m_column_names && m_column_types == other.m_column_types ;
	}

	bool CohortIndividualSource::ColumnSpec::operator!=( ColumnSpec const& other ) {
		return m_column_names != other.m_column_names || m_column_types != other.m_column_types ;
	}
	
	void CohortIndividualSource::ColumnSpec::add_column( std::string const& name, ColumnType const type ) {
		m_column_names.push_back( name ) ;
		m_column_types.push_back( type ) ;
	}
	
	CohortIndividualSource::ColumnSpec CohortIndividualSource::ColumnSpec::operator+( ColumnSpec const& other ) {
		std::vector< std::string > column_names = m_column_names ;
		std::vector< ColumnType > column_types = m_column_types ;
		column_names.insert( column_names.end(), other.m_column_names.begin(), other.m_column_names.end() ) ;
		column_types.insert( column_types.end(), other.m_column_types.begin(), other.m_column_types.end() ) ;
		return ColumnSpec( column_names, column_types ) ;
	}
	
	void CohortIndividualSource::ColumnSpec::remove( std::string const& name ) {
		std::size_t i = find_column( name ) ;
		m_column_names.erase( m_column_names.begin() + i ) ;
		m_column_types.erase( m_column_types.begin() + i ) ;
	}
	
	CohortIndividualSource::UniquePtr CohortIndividualSource::create(
		std::string source_spec,
		std::string const& missing_value,
		std::string const& choice
	) {
		std::size_t pos = source_spec.find( "://" ) ;
		if( pos == std::string::npos ) {
			source_spec = "file://" + source_spec ;
		}

		if( source_spec.substr( 0, 7 ) == "data://" ) {
			std::istringstream istr( source_spec.substr( 7, source_spec.size() )) ;
			if( choice == "strict" ) {
				return UniquePtr( new TraditionalStrictCohortIndividualSource( istr, missing_value )) ;
			}
			else if( choice == "categorical" ) {
				return UniquePtr( new CategoricalCohortIndividualSource( istr, missing_value )) ;
			}
		}
		else if( source_spec.substr( 0, 7 ) == "file://" ){
			std::string const filename = source_spec.substr( 7, source_spec.size() ) ;
			if( choice == "strict" ) {
				return UniquePtr( new TraditionalStrictCohortIndividualSource( filename, missing_value )) ;
			}
			else if( choice == "categorical" ) {
				return UniquePtr( new CategoricalCohortIndividualSource( filename, missing_value )) ;
			}
			else {
				assert(0) ;
			}
		}
		assert(0) ;
	}
	
	CohortIndividualSource const& CohortIndividualSource::get_base_source() const {
		return *this ;
	}

	CohortIndividualSource const& CohortIndividualSource::get_parent_source() const {
		return *this ;
	}
	
	std::string CohortIndividualSource::get_source_spec() const {
		return "(unknown)" ;
	}
	
	std::vector< std::size_t > CohortIndividualSource::find_entries( Entry const& entry, std::string const& column_name ) const {
		std::vector< std::size_t > result ;
		for( std::size_t i = 0; i < get_number_of_individuals(); ++i ) {
			if( get_entry( i, column_name ) == entry ) {
				result.push_back( i ) ;
			}
		}
		return result ;
	}
	
	void CohortIndividualSource::get_column_values( std::string const& column_name, boost::function< void ( std::size_t, VariantEntry ) > callback ) const {
		for( std::size_t i = 0; i < get_number_of_individuals(); ++i ) {
			callback( i, get_entry( i, column_name )) ;
		}
	}
	
	std::size_t CohortIndividualSource::ColumnSpec::get_number_of_covariates() const {
		return std::count( m_column_types.begin(), m_column_types.end(), e_DISCRETE_COVARIATE )
			+ std::count( m_column_types.begin(), m_column_types.end(), e_CONTINUOUS_COVARIATE ) ;
	}
	std::size_t CohortIndividualSource::ColumnSpec::get_number_of_phenotypes() const {
		return std::count( m_column_types.begin(), m_column_types.end(), e_BINARY_PHENOTYPE )
			+ std::count( m_column_types.begin(), m_column_types.end(), e_CONTINUOUS_PHENOTYPE ) ;
	}
	
	bool CohortIndividualSource::check_for_column( std::string const& column_name ) const {
		return get_column_spec().check_for_column( column_name ) ;
	}
	
	std::size_t CohortIndividualSource::get_number_of_covariates() const {
		return get_column_spec().get_number_of_covariates() ;
	}
	std::size_t CohortIndividualSource::get_number_of_phenotypes() const {
		return get_column_spec().get_number_of_phenotypes() ;
	}
	
	bool CohortIndividualSource::ColumnSpec::check_for_column( std::string const& column_name ) const {
		return std::find( m_column_names.begin(), m_column_names.end(), column_name ) != m_column_names.end() ;
	}

	std::size_t CohortIndividualSource::ColumnSpec::find_column( std::string const& column_name ) const {
		std::vector< std::string >::const_iterator
			where = std::find( m_column_names.begin(), m_column_names.end(), column_name ) ;
		if( where == m_column_names.end() ) {
			throw BadArgumentError( "CohortIndividualSource::ColumnSpec::find_column()", "column_name = \"" + column_name + "\".\n" ) ;
		}
		return std::size_t( where - m_column_names.begin() ) ;
	}
	
	std::ostream& operator<< ( std::ostream& out, CohortIndividualSource::ColumnType const& type ) {
		switch( type ) {
			case CohortIndividualSource::e_ID_COLUMN:
			case CohortIndividualSource::e_MISSINGNESS_COLUMN:
				out << "0" ;
				break ;
			case CohortIndividualSource::e_DISCRETE_COVARIATE:
				out << "D" ;
				break ;
			case CohortIndividualSource::e_CONTINUOUS_COVARIATE:
				out << "C" ;
				break ;
			case CohortIndividualSource::e_BINARY_PHENOTYPE:
				out << "B" ;
				break ;
			case CohortIndividualSource::e_CONTINUOUS_PHENOTYPE:
				out << "P" ;
				break ;
		}
		return out ;
	}

	std::ostream& operator<<( std::ostream& ostr, CohortIndividualSource::SingleColumnSpec const& spec ) {
		return ostr << spec.name() << ":" << spec.type() ;
	}

	std::ostream& operator<<( std::ostream& ostr, CohortIndividualSource::ColumnSpec const& spec ) {
		for( std::size_t i = 0; i < spec.size(); ++i ) {
			if( i > 0 ) {
				ostr << "," ;
			}
			ostr << spec[i] ;
		}
		return ostr ;
	}
	
}
