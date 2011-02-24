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
	
	CohortIndividualSource::ColumnSpec CohortIndividualSource::ColumnSpec::operator+( ColumnSpec const& other ) {
		std::vector< std::string > column_names = m_column_names ;
		std::vector< ColumnType > column_types = m_column_types ;
		column_names.insert( column_names.end(), other.m_column_names.begin(), other.m_column_names.end() ) ;
		column_types.insert( column_types.end(), other.m_column_types.begin(), other.m_column_types.end() ) ;
		return ColumnSpec( column_names, column_types ) ;
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
	
	
	std::size_t CohortIndividualSource::ColumnSpec::get_number_of_covariates() const {
		return std::count( m_column_types.begin(), m_column_types.end(), e_DISCRETE_COVARIATE )
			+ std::count( m_column_types.begin(), m_column_types.end(), e_CONTINUOUS_COVARIATE ) ;
	}
	std::size_t CohortIndividualSource::ColumnSpec::get_number_of_phenotypes() const {
		return std::count( m_column_types.begin(), m_column_types.end(), e_BINARY_PHENOTYPE )
			+ std::count( m_column_types.begin(), m_column_types.end(), e_CONTINUOUS_PHENOTYPE ) ;
	}
	
	bool CohortIndividualSource::check_for_column( std::string const& column_name ) const {
		std::vector< std::string > column_names = get_column_spec().get_names() ;
		return std::find( column_names.begin(), column_names.end(), column_name ) != column_names.end() ;
	}
	
	std::size_t CohortIndividualSource::get_number_of_covariates() const {
		return get_column_spec().get_number_of_covariates() ;
	}
	std::size_t CohortIndividualSource::get_number_of_phenotypes() const {
		return get_column_spec().get_number_of_phenotypes() ;
	}
	
	
	CohortIndividualSource::Entry::Entry():
		m_entrydata( MissingValue() )
	{}
	
	bool CohortIndividualSource::Entry::is_missing() const {
		return boost::get< MissingValue >( &m_entrydata ) ;
	}
	
	template<> double CohortIndividualSource::Entry::as() const {
		if( int const* v = boost::get< int >( &m_entrydata )) {
			return double( *v ) ;
		}
		else {
			return boost::get< double >( m_entrydata ) ;
		} 
	}
	
	bool CohortIndividualSource::Entry::operator==( Entry const& rhs ) const {
		return m_entrydata == rhs.m_entrydata ;
	}
	
	bool CohortIndividualSource::Entry::operator<( Entry const& rhs ) const {
		return m_entrydata < rhs.m_entrydata ;
	}
	
	std::ostream& operator<<( std::ostream& ostr, CohortIndividualSource::Entry const& entry ) {
		return ostr << entry.m_entrydata ;
	}
}
