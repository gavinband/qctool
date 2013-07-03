
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <boost/regex.hpp>
#include "FlatFileFrequentistGenomeWideAssociationResults.hpp"

namespace {
	double const NA = std::numeric_limits< double >::quiet_NaN() ;
}

/*
* Class FlatFileFrequentistGenomeWideAssociationResults.
* This is a utility base class for classes reading scan results flat files.
* It provides:
* - storage for values read.
* - 
*
* - and implements the FrequentistGenomeWideAssociationResults interface for retrieving tresults..
*/
void FlatFileFrequentistGenomeWideAssociationResults::add_data(
	std::vector< genfile::wildcard::FilenameMatch > const& filenames,
	SNPResultCallback callback,
	ProgressCallback progress_callback
) {
	setup( filenames, callback, progress_callback ) ;
}

std::size_t FlatFileFrequentistGenomeWideAssociationResults::get_number_of_SNPs() const {
	return m_snps.size() ;
}

genfile::SNPIdentifyingData2 const& FlatFileFrequentistGenomeWideAssociationResults::get_SNP( std::size_t snp_i ) const {
	assert( snp_i < m_snps.size() ) ;
	return m_snps[ snp_i ] ;
}
void FlatFileFrequentistGenomeWideAssociationResults::get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const {
	*result = m_betas.row( snp_i ).cast< double >() ;
}
void FlatFileFrequentistGenomeWideAssociationResults::get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const {
	*result = m_ses.row( snp_i ).cast< double >() ;
}
void FlatFileFrequentistGenomeWideAssociationResults::get_pvalue( std::size_t snp_i, double* result ) const {
	*result = m_pvalues( snp_i ) ;
}
void FlatFileFrequentistGenomeWideAssociationResults::get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const {
	*result = m_sample_counts.row( snp_i ).cast< double >() ;
}
void FlatFileFrequentistGenomeWideAssociationResults::get_info( std::size_t snp_i, double* result ) const {
	*result = m_info( snp_i ) ;
}
void FlatFileFrequentistGenomeWideAssociationResults::get_maf( std::size_t snp_i, double* result ) const {
	*result = m_maf( snp_i ) ;
}
void FlatFileFrequentistGenomeWideAssociationResults::get_frequency( std::size_t snp_i, double* result ) const {
	*result = ( 2.0 * m_sample_counts( snp_i, 2 ) + m_sample_counts( snp_i, 1 ) ) / ( 2.0 * m_sample_counts.row( snp_i ).sum() ) ;
}
void FlatFileFrequentistGenomeWideAssociationResults::get_variable( std::size_t snp_i, std::string const& variable, double* value ) const {
	std::map< std::string, std::vector< double > >::const_iterator where = m_extra_variables.find( variable ) ;
	assert( where != m_extra_variables.end() ) ;
	*value = where->second[ snp_i ] ;
}

std::string FlatFileFrequentistGenomeWideAssociationResults::get_summary( std::string const& prefix, std::size_t target_column ) const {
	using genfile::string_utils::to_string ;
	
	// estimate memory used in SNPs.
	unsigned long mem_used = 0 ;
	for( std::size_t i = 0; i < m_snps.size(); ++i ) {
		mem_used += m_snps[i].get_estimated_bytes_used() ;
	}
	
	std::string result = "scan results object ("
		+ to_string( m_snps.size() )
		+ " SNPs"
		+ ", ~"
		+ to_string(
			(
				(
					m_betas.size()
					+ m_ses.size()
					+ m_pvalues.size()
					+ m_info.size()
					+ m_maf.size()
					+ m_sample_counts.size()
				) * sizeof( float )
				+ mem_used
			)
			/ 1000000.0
		) + "Mb in use.)" ;
	return result ;
}

FlatFileFrequentistGenomeWideAssociationResults::FlatFileFrequentistGenomeWideAssociationResults():
	m_missing_value( "NA" )
{}

void FlatFileFrequentistGenomeWideAssociationResults::setup(
	std::vector< genfile::wildcard::FilenameMatch > const& filenames,
	SNPResultCallback callback,
	ProgressCallback progress_callback
) {
	progress_callback( 0, 100 ) ;
	setup( statfile::BuiltInTypeStatSource::open( filenames ), callback, progress_callback ) ;
}

void FlatFileFrequentistGenomeWideAssociationResults::setup(
	statfile::BuiltInTypeStatSource::UniquePtr source,
	SNPResultCallback callback,
	ProgressCallback progress_callback
) {
	ColumnMap column_map = get_columns_to_store( *source ) ;
	
	int degrees_of_freedom = ( column_map.left.find( "_beta_2" ) == column_map.left.end() ) ? 1 : 2 ;
	
	std::size_t snp_index = m_snps.size() ;

	genfile::SNPIdentifyingData snp ;
	Eigen::VectorXd betas ;
	Eigen::VectorXd ses ;
	
	for( ; read_snp( *source, snp ); (*source) >> statfile::ignore_all() ) {
		
		if( snp_index >= m_snps.size() ) {
			resize_storage( snp_index + 100000, degrees_of_freedom ) ;
		}
		
		// Deal with strange non-ids.  This isn't a general solution but 
		if(
			snp.get_SNPID() == "---" // IMPUTE2
			|| snp.get_SNPID() == "?" // QCTOOL under some usages
			|| snp.get_SNPID() == "NA" // dunno.
		) {
			snp.set_SNPID( "" ) ;
		}

		// read data columns

		ColumnMap::right_const_iterator i = column_map.right.begin(), end_i = column_map.right.end() ;
		std::string value ;
		for( ; i != end_i; ++i ) {
			(*source)
				>> statfile::ignore( i->first - source->current_column() )
				>> value ;

			store_value(
				snp_index,
				i->second,
				( value == m_missing_value ) ? NA : genfile::string_utils::to_repr< double >( value )
			) ;
		}

		m_snps[ snp_index ] = snp ;
		if( check_if_snp_accepted( snp_index ) ) {
			++snp_index ;
		}

		if( progress_callback ) {
			progress_callback( double( source->number_of_rows_read() + 1 ), source->number_of_rows() ) ;
		}
	}
	
	// Now deallocate any unused memory.
	m_snps.resize( snp_index ) ;
	resize_storage( snp_index, degrees_of_freedom ) ;
}

FlatFileFrequentistGenomeWideAssociationResults::ColumnMap FlatFileFrequentistGenomeWideAssociationResults::get_columns_to_store(
	statfile::BuiltInTypeStatSource const& source
) {
	using genfile::string_utils::to_string ;
	std::set< std::string > desired_columns = get_desired_columns() ;

	ColumnMap result ;
	for( std::size_t i = 0; i < source.number_of_columns(); ++i ) {
		std::string name = source.name_of_column( i ) ;
		for( std::set< std::string >::iterator j = desired_columns.begin(); j != desired_columns.end(); ++j ) {
			assert( j->size() > 0 ) ;
			boost::regex regex( *j ) ;
			if( boost::regex_match( name, regex ) ) {
				std::pair< ColumnMap::iterator, bool > insertion = result.insert( ColumnMap::value_type( *j, i )) ;
				if( !insertion.second ) {
					throw genfile::OperationFailedError(
						"SNPTESTResults::get_columns_to_store()",
						"m_column_map",
						"Insertion of value for column \"" + name + "\" (matching \"" + to_string( *j ) + "\")."
					) ;
				}
			}
		}
	}
	
	for( std::set< std::string >::const_iterator i = desired_columns.begin(); i != desired_columns.end(); ++i ) {
		if( result.left.find( *i ) == result.left.end() ) {
			throw genfile::BadArgumentError(
				"FlatFileFrequentistGenomeWideAssociationResults::get_columns_to_store()",
				"required column=\"" + *i + "\"",
				"Could not find column matching \"" + *i + "\" in source \"" + source.get_source_spec() + "\""
			) ;
		}
	}

	return result ;
}

void FlatFileFrequentistGenomeWideAssociationResults::resize_storage( Eigen::MatrixXf::Index const N_snps, Eigen::MatrixXf::Index const degrees_of_freedom ) {
	using std::min ;
	int const current_N = min( N_snps, Eigen::MatrixXf::Index( m_snps.size() ) ) ;
	{
		// free any unused memory for SNPs.
		std::vector< genfile::SNPIdentifyingData2 > snps( N_snps ) ;
		std::copy( m_snps.begin(), m_snps.end(), snps.begin() ) ;
		m_snps.swap( snps ) ;
	}
	{
		Eigen::MatrixXf betas = Eigen::MatrixXf::Zero( N_snps, degrees_of_freedom ) ;
		if( m_betas.rows() > 0 ) {
			betas.block( 0, 0, current_N, degrees_of_freedom ) = m_betas.block( 0, 0, current_N, degrees_of_freedom ) ;
		}
		m_betas.swap( betas ) ;
	}
	{
		Eigen::MatrixXf ses = Eigen::MatrixXf::Zero( N_snps, degrees_of_freedom )  ;
		if( m_ses.rows() > 0 ) {
			ses.block( 0, 0, current_N, degrees_of_freedom ) = m_ses.block( 0, 0, current_N, degrees_of_freedom ) ;
		}
		m_ses.swap( ses ) ;
	}
	{
		Eigen::VectorXf pvalues = Eigen::VectorXf::Zero( N_snps ) ;
		pvalues.head( current_N ) = m_pvalues.head( current_N ) ;
		m_pvalues.swap( pvalues ) ;
	}
	{
		Eigen::VectorXf info = Eigen::VectorXf::Zero( N_snps ) ;
		info.head( current_N ) = m_info.head( current_N ) ;
		m_info.swap( info ) ;
	}
	{
		Eigen::VectorXf maf = Eigen::VectorXf::Zero( N_snps ) ;
		maf.head( current_N ) = m_maf.head( current_N ) ;
		m_maf.swap( maf ) ;
	}
	{
		Eigen::MatrixXf sample_counts = Eigen::MatrixXf::Zero( N_snps, 4 ) ;
		if( m_sample_counts.rows() > 0 ) {
			sample_counts.block( 0, 0, current_N, 4 ) = m_sample_counts.block( 0, 0, current_N, 4 ) ;
		}
		m_sample_counts.swap( sample_counts ) ;
	}
	{
		for( ExtraVariables::iterator i = m_extra_variables.begin(); i != m_extra_variables.end(); ++i ) {
			i->second.resize( N_snps ) ; // std::vector resize does not lose data.
			std::vector< double > v = i->second ;
			i->second.swap( v ) ;
		}
	}
}

void FlatFileFrequentistGenomeWideAssociationResults::free_unused_memory() {
	std::size_t const N_snps = m_snps.size() ;
	{
		std::vector< genfile::SNPIdentifyingData2 > snps( m_snps.begin(), m_snps.end() ) ;
		m_snps.swap( snps ) ;
	}
	{
		Eigen::MatrixXf betas = m_betas.block( 0, 0, N_snps, m_betas.cols() ) ;
		m_betas.swap( betas ) ;
	}
	{
		Eigen::MatrixXf ses = m_ses.block( 0, 0, N_snps, m_betas.cols() ) ; ;
		m_ses.swap( ses ) ;
	}
	{
		Eigen::VectorXf pvalues = m_pvalues.head( N_snps ) ;
		m_pvalues.swap( pvalues ) ;
	}
	{
		Eigen::VectorXf info = m_info.head( N_snps ) ;
		m_info.swap( info ) ;
	}
	{
		Eigen::VectorXf maf = m_maf.head( N_snps ) ;
		m_maf.swap( maf ) ;
	}
	{
		Eigen::MatrixXf sample_counts = m_sample_counts.block( 0, 0, N_snps, m_sample_counts.cols() ) ;
		m_sample_counts.swap( sample_counts ) ;
	}
	{
		for( ExtraVariables::iterator i = m_extra_variables.begin(); i != m_extra_variables.end(); ++i ) {
			i->second.resize( N_snps ) ;
			std::vector< double > v = i->second ;
			i->second.swap( v ) ;
		}
	}
}

