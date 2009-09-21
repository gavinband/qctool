#include <set>
#include <cassert>
#include <fstream>
#include <iostream>

#include "GenRow.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SNPInListCondition.hpp"
#include "GenRowStatistics.hpp"
#include "string_utils.hpp"

#include "SNPDataSource.hpp"
#include "wildcard.hpp"

SNPInListCondition::SNPInListCondition( std::string const& filename )
 : m_filenames( std::size_t(1), filename )
{
	setup() ;
}

SNPInListCondition::SNPInListCondition( std::vector< std::string > const& filenames )
 : m_filenames( filenames )
{
	setup() ;
}

void SNPInListCondition::setup() {
	for( std::size_t i = 0; i < m_filenames.size(); ++i ) {
		if( file_appears_to_be_plain( m_filenames[i] )) {
			load_from_plain_file( m_filenames[i] ) ;			
		}
		else {
			load_from_gen_files( m_filenames[i] ) ;
		}
	}
}

bool SNPInListCondition::file_appears_to_be_plain( std::string const& filename ) {
	// return true if the first (up to) ten lines (or 1000-character blocks) contain a single
	// string containing only printable non-whitespace chars.
	std::ifstream file( filename.c_str() ) ;
	std::vector< char > buffer( 1000 ) ;

	for( std::size_t count = 0; file.get( &buffer[0], buffer.size(), '\n' ) && count < 10 ; ++count ) {
		std::size_t count = file.gcount() ;
		if( buffer[ count - 1] == '\n' ) {
			--count ;
		}
		std::string line( buffer.begin(), buffer.begin() + count ) ;
		std::size_t pos = line.find_first_not_of( "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-*:$" ) ;
		if( pos != std::string::npos ) {
			return false ;
		}
	}
	return true ;
}


void SNPInListCondition::load_from_gen_files( std::string const& filename ) {
	std::auto_ptr< genfile::SNPDataSource > source = genfile::SNPDataSource::create( wildcard::find_files_matching_path_with_integer_wildcard( filename )) ;
	uint32_t number_of_samples, SNP_position ;
	std::set< std::string > SNP_positions ;
	
	while( (*source).get_snp_identifying_data( 
			genfile::set_value( number_of_samples ),
			genfile::ignore(),
			genfile::ignore(),
			genfile::set_value( SNP_position ),
			genfile::ignore(),
			genfile::ignore()
		)
	) {
		SNP_positions.insert( to_string( SNP_position )) ;
		(*source).ignore_snp_probability_data( number_of_samples ) ;
	}

	if( source->number_of_snps_read() != source->total_number_of_snps() ) {
		throw genfile::FileStructureInvalidError() ;
	}
	m_id_list.insert( SNP_positions.begin(), SNP_positions.end() ) ;
}

void SNPInListCondition::load_from_plain_file( std::string const& filename ) {
	FromFileSet< std::set< std::string > > strings_in_file( filename ) ;
	m_id_list.insert( strings_in_file.begin(), strings_in_file.end() ) ;
}

bool SNPInListCondition::check_if_satisfied( string_to_value_map const& statistics ) const {
	GenRowStatistics const* row_statistics_ptr = dynamic_cast< GenRowStatistics const* >( &statistics ) ;
	if( !row_statistics_ptr ) {
		throw ConditionException( "SNPInListCondition only supports GenRowStatistics." ) ;
	}
	return list_contains( row_statistics_ptr->row().SNPID() )
		|| list_contains( row_statistics_ptr->row().RSID() )
		|| list_contains( to_string( row_statistics_ptr->row().SNP_position() ) ) ;
}

bool SNPInListCondition::list_contains( std::string const& elt ) const {
	return m_id_list.find( elt ) != m_id_list.end() ;
}

void SNPInListCondition::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< "snp-in-list(" ;
	for( std::size_t i = 0; i < m_filenames.size(); ++i ) {
		if( i > 0 ) {
			oStream << ", " ;
		}
		oStream << m_filenames[i] ;
	}
	oStream << ")" ;
}