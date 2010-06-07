#ifndef GENFILE_TRADITIONALSTRICTCOHORTINDIVIDUALSOURCE_HPP
#define GENFILE_TRADITIONALSTRICTCOHORTINDIVIDUALSOURCE_HPP

#include <vector>
#include <string>
#include <iosfwd>
#include "genfile/FromFileCohortIndividualSource.hpp"

namespace genfile {
	// class TraditionalCohortIndividualSource
	// This class implements a strict version of the sample file format described here:
	// http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format_new.html
	// Each of the id_1 and id_2 fields must contain unique entries.
	class TraditionalStrictCohortIndividualSource: public FromFileCohortIndividualSource
	{
	public:
		TraditionalStrictCohortIndividualSource( std::string const& filename, std::string const& missing_values = "NA" ) ;
		TraditionalStrictCohortIndividualSource( std::istream& stream, std::string const& missing_values = "NA" ) ;
	private:
		static Entry get_entry_from_string( std::string const& entry_as_string, ColumnType column_type ) ;
		void check_sample_ids() const ;
	} ;
}

#endif
