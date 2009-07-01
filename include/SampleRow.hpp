#ifndef __GTOOL__SAMPLEROW_HPP__
#define __GTOOL__SAMPLEROW_HPP__

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "GToolException.hpp"

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

class SampleRow
{
	public:

		// Read column headings and types from strings, in which the entities should be whitespace-separated.
		SampleRow() ;
		SampleRow( std::vector<std::string> const& column_headings, std::vector<char> const& column_types ) ;

	public:

		void reset( std::vector<std::string> const& column_headings, std::vector<char> const& column_types ) ;

		std::string ID1() const { return m_id1 ; }
		std::string ID2() const { return m_id2 ; }

		double further_data( std::string const& ) const ;
		std::vector<std::string> const& column_headings() const { return m_column_headings ; }
		std::vector<char> const& column_types() const { return m_column_types ; }
		
		bool have_column( std::string const& heading ) const ;
		void add_column( std::string const& heading, char type ) ;
		void set_value( std::string const& heading, double value ) ;
		
		friend std::istream& operator>>( std::istream& aStream, SampleRow& row ) ;

	private:
		std::string m_id1, m_id2 ;
		std::map< std::string, double > m_further_data ;
		std::vector<std::string> m_column_headings ;
		std::vector<char> m_column_types ;		

		SampleRow( SampleRow const& other ) ;
} ;

std::ostream& operator<<( std::ostream& aStream, SampleRow const& row ) ;

#endif
