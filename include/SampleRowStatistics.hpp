
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef __GTOOL_SAMPLEROWSTATISTICS__
#define __GTOOL_SAMPLEROWSTATISTICS__

#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "SampleRow.hpp"
#include "GenotypeProportions.hpp"
#include "GToolException.hpp"
#include "GenotypeAssayStatistics.hpp"

struct SampleRowStatistics: public GenotypeAssayStatistics
{
	typedef GenotypeAssayStatistics base_t ;
	public:
		SampleRowStatistics() {} ;
		void process( SampleRow const& row, GenotypeProportions const& amounts, std::size_t number_of_samples ) ;

		SampleRow const& row() const { return *m_row ; }

		void add_to_sample_row( SampleRow& row, std::string const&, std::string = "" ) const ;

	private:
		
		SampleRow const* m_row ;
} ;

struct SampleRowSpecificStatistic: public GenotypeAssayStatistic
{
	typedef GenotypeAssayStatistic base_t ;
	public:
		virtual double calculate_value( SampleRowStatistics const& row_statistics ) const = 0;
		virtual std::string calculate_string_value( SampleRowStatistics const& row_statistics ) const ;
	
		// These methods forward to the ones above.  An exception is thrown if the statistics are not
		// of type SampleRowStatistics.
		double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
		std::string calculate_string_value( GenotypeAssayStatistics const& statistics ) const ;

	protected:
		SampleRowStatistics const& get_row_statistics( GenotypeAssayStatistics const& statistics ) const ;
} ;

struct SampleRowID1: public SampleRowSpecificStatistic
{
	double calculate_value( SampleRowStatistics const& row_statistics ) const ;
	std::string calculate_string_value( SampleRowStatistics const& row_statistics ) const ;
} ;

struct SampleRowID2: public SampleRowSpecificStatistic
{
	double calculate_value( SampleRowStatistics const& row_statistics ) const ;
	std::string calculate_string_value( SampleRowStatistics const& row_statistics ) const ;
} ;

#endif

