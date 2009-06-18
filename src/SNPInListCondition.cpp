#include <set>
#include <cassert>
#include <fstream>

#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SNPInList.hpp"
#include "Whitespace.hpp"

SNPInList::SNPInList( std::string filename )
	: m_id_list( filename ) 
{}

bool SNPInList::check_if_satisfied( GenRow const& genRow, GenotypeAssayStatistics const* ) const {
	return m_id_list.find( genRow.SNPID() ) != m_id_list.end() ;
}
