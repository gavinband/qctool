#ifndef __GTOOL_GENOTYPEASSAYSTATISTICS__
#define __GTOOL_GENOTYPEASSAYSTATISTICS__

#include <vector>
#include <iostream>
#include <map>
#include "GenotypeProportions.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "GenotypeAssayBasicStatistics.hpp"

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

// This class holds data representing an assay of one or more SNPs.
// It has methods to return the amounts of AA, AB and BB genotypes in the sample,
// as well as genotype proportions and allele proportions. 
struct GenotypeAssayStatistics: public GenotypeAssayBasicStatistics
{
	typedef GenotypeAssayBasicStatistics base_t ;

	public:
		GenotypeAssayStatistics() ;
		~GenotypeAssayStatistics() ;

		template< typename Iterator >
		void process( Iterator begin, Iterator const& end ) {
			reset_statistics() ;
			base_t::process( begin, end ) ;
		}
		
		// Methods to manipulate list of statistics
		void add_statistic( std::string const& name, std::auto_ptr< GenotypeAssayStatistic > statistic_ptr ) ;
		template< typename T >
		T get_statistic_value( std::string const& name ) const ;

	private:

		void reset_statistics() ;
	
		typedef std::map< std::string, GenotypeAssayStatistic* > statistics_t ;
		statistics_t m_statistics ;

	public:
		std::ostream& format_column_headers( std::ostream& ) ;
		std::ostream& format_statistic_values( std::ostream& aStream ) const ;
} ;

// Base class for individual statistics
// This class contains a statistic value cache.
// The reset() method is provided to reset this cache before
// using the statistic on a different GenotypeAssayStatistics object.
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

