
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <boost/regex.hpp>
#include <boost/format.hpp>
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
void FlatFileFrequentistGenomeWideAssociationResults::get_covariance_upper_triangle( std::size_t snp_i, Eigen::VectorXd* result ) const {
	*result = m_covariance.row( snp_i ).cast< double >() ;
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
void FlatFileFrequentistGenomeWideAssociationResults::get_variable( std::size_t snp_i, std::string const& variable, std::string* value ) const {
	std::map< std::string, std::vector< std::string > >::const_iterator where = m_extra_variables.find( variable ) ;
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
	
	std::ostringstream str( result ) ;
	str << "\n" << " -- using the following columns:\n" ;
	SourceColumnMap::right_const_iterator
		i = m_column_map.get().right.begin(),
		end_i = m_column_map.get().right.end() ;
	
	for( std::size_t count = 0; i != end_i; ++i, ++count ) {
		str << " -- column " << i->first << ": " << "(\"" << i->second.name() << "\", matching \"" << i->second.pattern().str() << "\")\n" ;
	}
	return str.str() ;
}

FlatFileFrequentistGenomeWideAssociationResults::FlatFileFrequentistGenomeWideAssociationResults():
	m_missing_value( "NA" ),
	m_degrees_of_freedom( 0 )
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
	if( !m_column_map ) {
		m_desired_columns = setup_columns( source->column_names() ) ;
		m_column_map = get_source_column_map( *source, m_desired_columns ) ;
	} else {
		check_columns( *source, m_column_map.get() ) ;
	}
	assert( m_column_map ) ;
	int const degrees_of_freedom = get_effect_parameter_names().size() ;
	
	std::size_t snp_index = m_snps.size() ;

	genfile::SNPIdentifyingData snp ;
	
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

		// read data columns.
		// We use the index-to-regex map to iterate over required columns in the right order.
		SourceColumnMap::right_const_iterator i = m_column_map.get().right.begin(), end_i = m_column_map.get().right.end() ;
		std::string value ;
		for( ; i != end_i; ++i ) {
			(*source)
				>> statfile::ignore( i->first - source->current_column() )
				>> value ;

			store_value(
				snp_index,
				i->second.simplified_name(),
				value
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

FlatFileFrequentistGenomeWideAssociationResults::SourceColumnMap FlatFileFrequentistGenomeWideAssociationResults::get_source_column_map(
	statfile::BuiltInTypeStatSource const& source,
	DesiredColumns const& desired_columns
) const {
	SourceColumnMap result ;
	for( std::size_t i = 0; i < desired_columns.size(); ++i ) {
		std::pair< SourceColumnMap::iterator, bool > insertion = result.insert(
			SourceColumnMap::value_type(
				desired_columns[i], desired_columns[i].index()
			)
		) ;
		if( !insertion.second ) {
			throw genfile::OperationFailedError(
				"FlatFileFrequentistGenomeWideAssociationResults::get_source_column_map()",
				"m_column_map",
				"Insertion of value for column \"" + desired_columns[i].name() + "\" (matching \"" + desired_columns[i].pattern().str() + "\")."
			) ;
		}
	}
	return result ;
}

void FlatFileFrequentistGenomeWideAssociationResults::check_columns(
	statfile::BuiltInTypeStatSource const& source,
	SourceColumnMap const& column_map
) const {
	std::vector< std::string > const& column_names = source.column_names() ;
	// typedef boost::bimap< std::string, std::size_t > SourceColumnMap ;

	boost::format fmt( "Column %d (\"%s\") does not match expected name \"%s\"." ) ;
	SourceColumnMap::const_iterator
		i = column_map.begin(),
		end_i = column_map.end() ;
	for( ; i != end_i; ++i ) {
		std::string const& column_name = column_names[ i->right ] ;
		std::string const& expected_name = i->left.name() ;
		if( column_name != expected_name ) {
			throw genfile::BadArgumentError(
				"FlatFileFrequentistGenomeWideAssociationResults::check_columns()",
				"source=\"" + source.get_source_spec() + "\"",
				( fmt % ( i->right + 1 ) % column_name % expected_name ).str()
			) ;
		} ;
	}
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
		Eigen::MatrixXf covariance = Eigen::MatrixXf::Zero( N_snps, ( degrees_of_freedom - 1 ) * degrees_of_freedom / 2 )  ;
		if( m_covariance.rows() > 0 ) {
			covariance.block( 0, 0, current_N, covariance.cols() ) = m_covariance.block( 0, 0, current_N, covariance.cols() ) ;
		}
		m_covariance.swap( covariance ) ;
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
		Eigen::MatrixXf sample_counts = Eigen::MatrixXf::Zero( N_snps, 6 ) ;
		if( m_sample_counts.rows() > 0 ) {
			sample_counts.block( 0, 0, current_N, 6 ) = m_sample_counts.block( 0, 0, current_N, 6 ) ;
		}
		m_sample_counts.swap( sample_counts ) ;
	}
	{
		for( ExtraVariables::iterator i = m_extra_variables.begin(); i != m_extra_variables.end(); ++i ) {
			i->second.resize( N_snps ) ; // std::vector resize does not lose data.
			std::vector< std::string > v = i->second ;
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
		Eigen::MatrixXf covariance = m_covariance.block( 0, 0, N_snps, m_covariance.cols() ) ; ;
		m_covariance.swap( covariance ) ;
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
			std::vector< std::string > v = i->second ;
			i->second.swap( v ) ;
		}
	}
}

