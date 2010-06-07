#ifndef __GTOOL__SAMPLEROW_HPP__
#define __GTOOL__SAMPLEROW_HPP__

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "string_to_value_map.hpp"
#include "GToolException.hpp"
#include "genfile/CohortIndividualSource.hpp"

struct SampleRowException: public GToolException
{
	SampleRowException( std::string const& msg )
	: GToolException( msg )
	{}
} ;

struct BadSampleRowFormatException: public SampleRowException
{
	BadSampleRowFormatException( std::string const& msg )
	: SampleRowException( msg )
	{}
} ;


class SampleFileColumnTypes
{
		// 
		static char const eFirstEntries = '0' ;
		static char const eDiscreteCovariateForMantelHaentzelTest = '1' ;
		static char const eDiscreteCovariateForCombinedTest = '1' ;
		static char const eContinuousCovariate = '3' ;
		static char const ePhenotype = 'P' ;
} ;

class SampleRow: public string_to_value_map
{
	public:

		// Read column headings and types from strings, in which the entities should be whitespace-separated.
		SampleRow() ;
		//SampleRow( std::vector<std::string> const& column_headings, std::vector<char> const& column_types ) ;
		SampleRow( SampleRow const& other ) ;

		virtual ~SampleRow() {} ;
		
		typedef genfile::CohortIndividualSource::Entry Entry ;
	public:

		void reset( std::vector< genfile::CohortIndividualSource::SingleColumnSpec > const& column_spec ) ;

		std::string ID1() const ;
		std::string ID2() const ;

		Entry further_data( std::string const& ) const ;
		std::vector<std::string> const& column_headings() const { return m_column_headings ; }
		std::vector<char> const& column_types() const { return m_column_types ; }
		
		bool has_column( std::string const& heading ) const ;
		void add_column( std::string const& heading, char type, Entry value = 0.0 ) ;
		void set_value( std::string const& heading, Entry value ) ;

		//friend std::istream& operator>>( std::istream& aStream, SampleRow& row ) ;
		void read_ith_sample_from_source( std::size_t i, genfile::CohortIndividualSource const& source ) ;
	protected:
		
		// string_to_value_map derived methods.
		bool has_value( std::string const& name ) const ;
		double get_double_value( std::string const& name ) const ;
		std::string get_string_value( std::string const& name ) const ;

	private:
		typedef std::map< std::string, Entry > FurtherData ;
		FurtherData m_further_data ;
		std::vector<std::string> m_column_headings ;
		std::vector<char> m_column_types ;		

} ;

std::ostream& operator<<( std::ostream& aStream, SampleRow const& row ) ;

#endif
