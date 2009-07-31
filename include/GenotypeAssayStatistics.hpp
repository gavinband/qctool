#ifndef __GTOOL_GENOTYPEASSAYSTATISTICS__
#define __GTOOL_GENOTYPEASSAYSTATISTICS__

#include <vector>
#include <iostream>
#include <map>
#include "GenotypeProportions.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "GenotypeAssayBasicStatistics.hpp"
#include "string_to_value_map.hpp"

struct GenotypeAssayStatisticException: public GToolException
{
	GenotypeAssayStatisticException( std::string const& msg )
		: GToolException( msg )
	{}
};

struct StatisticNotFoundException: public GenotypeAssayStatisticException
{
	StatisticNotFoundException( std::string const& msg )
		: GenotypeAssayStatisticException( msg )
	{}
};


// Forward declaration
struct GenotypeAssayStatistic ;

struct GenotypeAssayStatistics: public GenotypeAssayBasicStatistics, public string_to_value_map
{
	typedef GenotypeAssayBasicStatistics base_t ;

	public:
		GenotypeAssayStatistics() ;
		~GenotypeAssayStatistics() ;

		// Methods to manipulate list of statistics
		void add_statistic( std::string const& name, std::auto_ptr< GenotypeAssayStatistic > statistic_ptr ) ;
		std::size_t size() const { return m_statistics.size() ; }

	protected:
		
		bool has_value( std::string const& name ) const ;
		double get_double_value( std::string const& name ) const ;
		std::string get_string_value( std::string const& name ) const ;
		
		typedef std::map< std::string, GenotypeAssayStatistic* > statistics_t ;
		typedef statistics_t::const_iterator statistic_iterator_t ;
		std::vector<std::string>::const_iterator begin_statistics() const { return m_statistic_names.begin() ; }
		std::vector<std::string>::const_iterator end_statistics() const { return m_statistic_names.end() ; }


	private:

		void reset() ;
	
		statistics_t m_statistics ;
		std::vector< std::string > m_statistic_names ;

	public:
		std::ostream& format_column_headers( std::ostream& ) ;
		std::ostream& format_statistic_values( std::ostream& aStream ) const ;
} ;

// Base class for individual statistics
struct GenotypeAssayStatistic
{
	public:
		GenotypeAssayStatistic() ;
		virtual ~GenotypeAssayStatistic() {}

		template< typename T>
		T get_value( GenotypeAssayStatistics const& ) const ;

		virtual void reset() const {};

		void set_precision( std::size_t precision ) { m_precision = precision ; }

	protected:
		virtual double calculate_value( GenotypeAssayStatistics const& ) const = 0;
		virtual std::string calculate_string_value( GenotypeAssayStatistics const& ) const ;
	
	private:
	
		std::size_t m_precision ;
} ;


#endif

