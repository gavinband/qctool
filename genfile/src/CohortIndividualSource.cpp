#include <vector>
#include <string>
#include <set>
#include <cassert>
#include <exception>
#include <iostream>
#include <fstream>
#include <boost/variant.hpp>

#include "genfile/Error.hpp"
#include "string_utils/string_utils.hpp"
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

	std::size_t CohortIndividualSource::get_number_of_columns() const {
		return 3 + get_number_of_covariates() + get_number_of_phenotypes() ;
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
	
	bool CohortIndividualSource::SingleColumnSpec::operator==( CohortIndividualSource::SingleColumnSpec const& right ) const {
		return Base( *this ) == Base( right ) ;
	}
	
	bool CohortIndividualSource::SingleColumnSpec::operator!=( CohortIndividualSource::SingleColumnSpec const& right ) const {
		return Base( *this ) != Base( right ) ;
	}
	
	CohortIndividualSource::UniquePtr CohortIndividualSource::create(
		std::string const& source_spec,
		std::string const& missing_value,
		std::string const& choice
	) {
		if( choice == "strict" ) {
			return UniquePtr( new TraditionalStrictCohortIndividualSource( source_spec, missing_value )) ;
		}
		else if( choice == "categorical" ) {
			return UniquePtr( new CategoricalCohortIndividualSource( source_spec, missing_value )) ;
		}
		else {
			assert(0) ;
		}
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
