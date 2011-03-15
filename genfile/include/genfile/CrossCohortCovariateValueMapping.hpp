#ifndef GENFILE_CROSS_COHORT_COVARIATE_MAPPING_HPP
#define GENFILE_CROSS_COHORT_COVARIATE_MAPPING_HPP

#include <vector>
#include <map>
#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	class CrossCohortCovariateValueMapping
		// Base classes for classes which map a given column of phenotypes or covariates,
		// across several cohorts, to transformed values.
	{
	public:
		typedef std::auto_ptr< CrossCohortCovariateValueMapping > UniquePtr ;
		static UniquePtr create( CohortIndividualSource::SingleColumnSpec const& column_spec ) ;

		typedef CohortIndividualSource::Entry Entry ;
	public:
		virtual ~CrossCohortCovariateValueMapping() {}
		virtual void add_source( CohortIndividualSource const& source ) = 0 ;
		virtual std::size_t get_number_of_unmapped_values() const = 0 ;
		virtual std::size_t get_number_of_distinct_unmapped_values() const = 0 ;
		virtual std::size_t get_number_of_distinct_mapped_values() const = 0 ;
		virtual std::size_t get_number_of_missing_values() const = 0 ;
		
		// Get the mapped value of the given entry.
		// Calling this with a value not included among the cohorts passed to this mapping
		// results in undefined behaviour.
		virtual Entry get_mapped_value( Entry const& level ) const = 0 ;

		// Get the unmapped value (the value occuring within the cohort sources) corresponding
		// to the given mapped value.  Calling this with a value not included among the cohorts
		// passed to this mapping results in undefined behaviour.
		virtual Entry get_unmapped_value( Entry const& level ) const = 0 ;

		// Return a human-readable summary of the mapping and/or its values.
		virtual std::string get_summary( std::string const& prefix = "" ) const = 0 ;
	} ;

	class LevelCountingCrossCohortCovariateValueMapping: public CrossCohortCovariateValueMapping
	{
	public:
		LevelCountingCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		void add_source( CohortIndividualSource const& source ) ;
		std::size_t get_number_of_unmapped_values() const ;
		std::size_t get_number_of_distinct_unmapped_values() const ;
		std::size_t get_number_of_missing_values() const ;
		
	public:
		typedef std::map< Entry, unsigned int > Entries ;

	protected:
		// function: entries()
		// Return a map from Entries to their multiplicities in the data.
		Entries const& entries() const { return m_entries ; }
		
	private:
		std::string m_column_name ;
		Entries m_entries ;
		std::size_t m_number_of_missing_values ;
		
	private:

		void add_entries_from_source( CohortIndividualSource const& source, std::string const& column_name ) ;
		// forbid copying / assignment
		LevelCountingCrossCohortCovariateValueMapping( LevelCountingCrossCohortCovariateValueMapping const& other ) ;
		LevelCountingCrossCohortCovariateValueMapping& operator=( LevelCountingCrossCohortCovariateValueMapping const& other ) ;
	} ;


	class CategoricalCrossCohortCovariateValueMapping: public LevelCountingCrossCohortCovariateValueMapping
		// A mapping which maps its values onto the inclusive range of positive integers
		// 1...number of distinct values.
	{
	public:
		CategoricalCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		std::size_t get_number_of_distinct_mapped_values() const ;

		// Return the entry corresponding to level i.
		// i must be a positive integer in the range 1...number of distinct entries
		Entry get_unmapped_value( Entry const& level ) const ;
		
		// Return the level (positive integer) corresponding to the given entry.
		// It is an error to pass this method an entry not present in the relevant column of the sources
		// this object has had added.
		Entry get_mapped_value( Entry const& entry ) const ;
		
		std::string get_summary( std::string const& prefix = "" ) const ;
		
	} ;
	
	class ContinuousVariableCrossCohortCovariateValueMapping: public LevelCountingCrossCohortCovariateValueMapping
	{
		// An identity mapping.  This is only useful insofar as it provides a nice summary of its data.
	public:
		ContinuousVariableCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		std::size_t get_number_of_distinct_mapped_values() const ;
		Entry get_unmapped_value( Entry const& level ) const ;
		Entry get_mapped_value( Entry const& entry ) const ;
		
		std::string get_summary( std::string const& prefix = "" ) const ;

	public:
		virtual std::string get_mapping_name() const { return "normalised" ; }

	protected:
		
		static void calculate_mean_and_variance( Entries const& entries, double* mean, double* variance ) ;
		typedef std::map< std::pair< double, double >, double > Histogram ;
		Histogram get_histogram(
			Entries const& entries,
		 	std::size_t number_of_bins
		) const ;
		std::string print_histogram(
			Histogram const& histogram,
			std::string const& prefix,
			std::size_t const height,
			std::size_t const bin_width
		) const ;
	} ;
	
	class NormalisingCrossCohortCovariateValueMapping: public ContinuousVariableCrossCohortCovariateValueMapping
		// A mapping which scales its given values (interpreted as real numbers) so that they have the given
		// mean and variance.
	{
	public:
		NormalisingCrossCohortCovariateValueMapping(
			std::string const& column_name
		) ;
		
		void add_source( CohortIndividualSource const& source ) ;
		
		std::size_t get_number_of_distinct_mapped_values() const ;
		Entry get_unmapped_value( Entry const& level ) const ;
		Entry get_mapped_value( Entry const& entry ) const ;
		
		// std::string get_summary( std::string const& prefix = "" ) const ;
		
		virtual std::string get_mapping_name() const { return "normalised" ; }
		
	private:
		double m_mean, m_variance, m_standard_deviation ;
	private:
		void analyse_values( Entries const& entries ) ;
	} ;
}

#endif
