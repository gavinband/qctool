#ifndef GENFILE_COHORTINDIVIDUALSOURCE_HPP
#define GENFILE_COHORTINDIVIDUALSOURCE_HPP

#include <string>
#include <memory>
#include <iosfwd>

#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>

#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/MissingValue.hpp"

namespace genfile {
	// Base class for classes which provide random-access view
	// to a set of samples, 
	class CohortIndividualSource
	{
	public:
		typedef std::auto_ptr< CohortIndividualSource > UniquePtr ;
		typedef std::auto_ptr< CohortIndividualSource const > ConstUniquePtr ;
		static UniquePtr create(
			std::string source_spec,
			std::string const& missing_value = "NA",
			std::string const& choice = "categorical"
		) ;
		struct ColumnSpec ;
	public:
		virtual ~CohortIndividualSource() {} ;
		virtual std::size_t get_number_of_individuals() const = 0 ;
		std::size_t get_number_of_covariates() const ;
		std::size_t get_number_of_phenotypes() const ;

		virtual ColumnSpec get_column_spec() const = 0 ;
		virtual bool check_for_column( std::string const& column_name ) const ;

		// method: get_entry()
		// get_entry returns the entry for the sample whose index is given and the named column.
		// is_missing(): return true if the value is missing, false otherwise
		// as< type >(): return the value.  type must either be int, double, or string, according to the type of the column.
		class Entry ;
		virtual Entry get_entry( std::size_t sample_i, std::string const& column_name ) const = 0 ;

		// Source objects may live in hierarchies.
		// This method return the eventual parent of the hierarchy.
		virtual CohortIndividualSource const& get_base_source() const ;
		// Find the parent of this source
		virtual CohortIndividualSource const& get_parent_source() const ;

		// method: get_source_spec()
		// get_source_spec() returns a human-readable specification for this source.
		virtual std::string get_source_spec() const ;

		// method find_entries()
		// find_entry returns the set of rows for which the given column equals the given entry.
		std::vector< std::size_t > find_entries( Entry const& entry, std::string const& column_name ) const ;

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
			bool is_phenotype() const ;
			bool is_covariate() const ;
			
			bool operator==( SingleColumnSpec const& right ) const ;
			bool operator!=( SingleColumnSpec const& right ) const ;
		} ;
		
		class ColumnSpec
		{
		public:
			ColumnSpec() ;
			ColumnSpec( std::vector< std::string > const& column_names, std::vector< ColumnType > const& column_types ) ;
			ColumnSpec( ColumnSpec const& other ) ;
			ColumnSpec& operator=( ColumnSpec const& other ) ;

			std::size_t size() const ;
			std::vector< std::string > get_names() const ;
			std::vector< ColumnType > get_types() const ;
			SingleColumnSpec get_spec( std::size_t i ) const ;
			SingleColumnSpec operator[]( std::size_t i ) const ;
			
			std::size_t get_number_of_covariates() const ;
			std::size_t get_number_of_phenotypes() const ;

			bool operator==( ColumnSpec const& other ) ;
			bool operator!=( ColumnSpec const& other ) ;
			
			// Concatenate two ColumnSpecs.
			ColumnSpec operator+( ColumnSpec const& other ) ;

		private:
			std::vector< std::string > m_column_names ;
			std::vector< ColumnType > m_column_types ;
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
