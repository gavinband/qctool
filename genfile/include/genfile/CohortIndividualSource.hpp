#ifndef GENFILE_COHORTINDIVIDUALSOURCE_HPP
#define GENFILE_COHORTINDIVIDUALSOURCE_HPP

#include <string>
#include <memory>
#include <iosfwd>

#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>

#include "genfile/Error.hpp"
#include "string_utils/string_utils.hpp"
#include "genfile/MissingValue.hpp"

namespace genfile {
	// Base class for classes which provide random-access view
	// to a set of samples, 
	class CohortIndividualSource
	{
	public:
		typedef std::auto_ptr< CohortIndividualSource > UniquePtr ;
		static UniquePtr create(
			std::string const& source_spec,
			std::string const& missing_value = "NA",
			std::string const& choice = "categorical"
		) ;
		struct SingleColumnSpec ;
	public:
		virtual ~CohortIndividualSource() {} ;
		virtual std::size_t get_number_of_columns() const ;
		virtual std::size_t get_number_of_individuals() const = 0 ;
		virtual std::size_t get_number_of_covariates() const = 0 ;
		virtual std::size_t get_number_of_phenotypes() const = 0 ;

		virtual std::vector< SingleColumnSpec > get_column_spec() const = 0 ;

		// method: get_entry()
		// get_entry returns the entry for the sample whose index is given and the named column.
		// is_missing(): return true if the value is missing, false otherwise
		// as< type >(): return the value.  type must either be int, double, or string, according to the type of the column.
		class Entry ;
		virtual Entry get_entry( std::size_t sample_i, std::string const& column_name ) const = 0 ;

	public:
		enum ColumnType { e_ID_COLUMN = 0, e_MISSINGNESS_COLUMN, e_DISCRETE_COVARIATE, e_CONTINUOUS_COVARIATE, e_BINARY_PHENOTYPE, e_CONTINUOUS_PHENOTYPE } ;
		friend std::ostream& operator<< ( std::ostream& out, ColumnType const& type ) ;

	public:
		
		class Entry {
		public:
			template< typename T > Entry( T const& value ) ;
			template< typename T > Entry& operator=( T const& value ) ;
			Entry() ; // Initialise with MissingValue.
		public:
			bool is_missing() const ;
			template< typename T > T as() const ;
		public:
			Entry( boost::variant< std::string, int, double, MissingValue > data ) ;

			bool operator==( Entry const& rhs ) const ;
			template< typename T > void operator==( T const& rhs ) const ; // this prevents implicit conversion to Entry.
			bool operator<( Entry const& rhs ) const ;
			template< typename T > void operator<( T const& rhs ) const ; // this prevents implicit conversion to Entry.

			friend std::ostream& operator<<( std::ostream&, Entry const& ) ;
		private:
			typedef boost::variant< MissingValue, std::string, int, double > EntryData ;
			EntryData m_entrydata ;
		} ;
		
		struct SingleColumnSpec: private std::pair< std::string, ColumnType >
		{
		private:
			typedef std::pair< std::string, ColumnType > Base ;
		public:
			SingleColumnSpec( std::string const& name = "", ColumnType const& type = e_ID_COLUMN ) ;
			SingleColumnSpec( SingleColumnSpec const& other ) ;
			SingleColumnSpec& operator=( SingleColumnSpec const& other ) ;

			std::string const& name() const ;
			ColumnType const type() const ;
			bool is_discrete() const ;
			bool is_continuous() const ;
			
			bool operator==( SingleColumnSpec const& right ) const ;
			bool operator!=( SingleColumnSpec const& right ) const ;
		} ;
	} ;

	template< typename T >
	CohortIndividualSource::Entry::Entry( T const& value ):
		m_entrydata( value )
	{}

	template< typename T >
	CohortIndividualSource::Entry& CohortIndividualSource::Entry::operator=( T const& value ) {
		this->m_entrydata = value ;
		return *this ;
	}
	
	template< typename T >
	T CohortIndividualSource::Entry::as() const {
		return boost::get< T >( m_entrydata ) ;
	}

	// The following specialisation allows the user to get either a double or integer entry as a double.
	// This simplifies some uses.
	template<> double CohortIndividualSource::Entry::as() const ;
}

#endif
