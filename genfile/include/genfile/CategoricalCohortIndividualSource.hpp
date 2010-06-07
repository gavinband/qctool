#ifndef GENFILE_CATEGORICALCOHORTINDIVIDUALSOURCE_HPP
#define GENFILE_CATEGORICALCOHORTINDIVIDUALSOURCE_HPP

#include <vector>
#include <string>
#include <iosfwd>
#include "genfile/FromFileCohortIndividualSource.hpp"

namespace genfile {
	// class CategoricalCohortIndividualSource
	// This class implements a version of the sample file format described here:
	// http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format_new.html
	// However, this class allows columns of type "D" (discrete) to contain entries other than
	// positive integers.  These columns are interpreted as strings.  The intention is to use this 
	// with an additional mechanism, such as the CrossCohortCovariateValueMapping class,
	// to map these values to categories.
	//
	// In addition, for Binary phenotype columns, we allow the coding "case/control" as well as "0/1"
	class CategoricalCohortIndividualSource: public FromFileCohortIndividualSource
	{
	public:
		CategoricalCohortIndividualSource( std::string const& filename, std::string const& missing_values = "NA" ) ;
		CategoricalCohortIndividualSource( std::istream& stream, std::string const& missing_values = "NA" ) ;
	private:
		static Entry get_entry_from_string( std::string const& entry_as_string, ColumnType column_type ) ;
		void check_sample_ids() const ;
	} ;
}

#endif
