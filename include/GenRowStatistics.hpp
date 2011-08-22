#ifndef __GTOOL_GENROWSTATISTICS__
#define __GTOOL_GENROWSTATISTICS__

#include <vector>
#include <iostream>
#include <map>
#include "GenRow.hpp"
#include "GenotypeProportions.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "GenotypeAssayStatistics.hpp"

struct GenRowStatistics: public GenotypeAssayStatistics
{
	typedef GenotypeAssayStatistics base_t ;
	public:
		GenRowStatistics() ;
		void process( GenRow const& row ) ;

		GenRow const& row() const { return *m_row ; }
	private:
	
		GenRow const* m_row ;
} ;

struct GenRowSpecificStatistic: public GenotypeAssayStatistic
{
	typedef GenotypeAssayStatistic base_t ;
	public:
		virtual double calculate_value( GenRowStatistics const& row_statistics ) const = 0;
		virtual std::string calculate_string_value( GenRowStatistics const& row_statistics ) const ;
	
		// These methods forward to the ones above.  An exception is thrown if the statistics are not
		// of type GenRowStatistics.
		double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
		std::string calculate_string_value( GenotypeAssayStatistics const& statistics ) const ;

	protected:
		GenRowStatistics const& get_row_statistics( GenotypeAssayStatistics const& statistics ) const ;
} ;


// Return Chromosome
struct GenRowChromosome: public GenRowSpecificStatistic
{
	double calculate_value( GenRowStatistics const& row_statistics ) const ;
	std::string calculate_string_value( GenRowStatistics const& row_statistics ) const ;
} ;

// Return SNP Position
struct GenRowSNPPosition: public GenRowSpecificStatistic
{
	double calculate_value( GenRowStatistics const& row_statistics ) const ;
	std::string calculate_string_value( GenRowStatistics const& row_statistics ) const ;
} ;

// Return SNPID
struct GenRowSNPID: public GenRowSpecificStatistic
{
	double calculate_value( GenRowStatistics const& row_statistics ) const ;
	std::string calculate_string_value( GenRowStatistics const& row_statistics ) const ;
} ;

// Return RSID
struct GenRowRSID: public GenRowSpecificStatistic
{
	double calculate_value( GenRowStatistics const& row_statistics ) const ;
	std::string calculate_string_value( GenRowStatistics const& row_statistics ) const ;
} ;

// Return minor or major allele
struct GenRowAllele: public GenRowSpecificStatistic
{
public:
	enum AlleleSelector { minor_allele = 0, major_allele = 1, first_allele = 2, second_allele = 3 } ;
	GenRowAllele( AlleleSelector const& selector ) ;
	double calculate_value( GenRowStatistics const& row_statistics ) const ;
	std::string calculate_string_value( GenRowStatistics const& row_statistics ) const ;
private:	
	AlleleSelector m_selector ;
} ;

#endif

