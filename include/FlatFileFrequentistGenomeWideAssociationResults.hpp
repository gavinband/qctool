
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef FLAT_FILE_FREQUENTISTGENOMEWIDEASSOCIATIONRESULTS_HPP
#define FLAT_FILE_FREQUENTISTGENOMEWIDEASSOCIATIONRESULTS_HPP

#include <boost/bimap.hpp>
#include "statfile/BuiltInTypeStatSource.hpp"
#include "FrequentistGenomeWideAssociationResults.hpp"

/*
* Class FlatFileFrequentistGenomeWideAssociationResults.
* This is a utility base class for classes reading scan results flat files.
* It provides:
* - storage for values read.
* - 
*
* - and implements the FrequentistGenomeWideAssociationResults interface for retrieving tresults..
*/
struct FlatFileFrequentistGenomeWideAssociationResults:  public FrequentistGenomeWideAssociationResults {
public:
	typedef std::auto_ptr< FlatFileFrequentistGenomeWideAssociationResults > UniquePtr ;
public:
	virtual void set_effect_size_column_regex( std::string const& beta_column_regex ) = 0 ;

	void add_data(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) ;

	std::size_t get_number_of_SNPs() const ;

	genfile::SNPIdentifyingData2 const& get_SNP( std::size_t snp_i ) const ;
	virtual int get_number_of_effect_parameters() const = 0 ;
	void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const ;
	void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const ;
	void get_covariance_upper_triangle( std::size_t snp_i, Eigen::VectorXd* result ) const ; 
	void get_pvalue( std::size_t snp_i, double* result ) const ;
	void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const ;
	void get_info( std::size_t snp_i, double* result ) const ;
	void get_maf( std::size_t snp_i, double* result ) const ;
	void get_frequency( std::size_t snp_i, double* result ) const ;
	void get_variable( std::size_t snp_i, std::string const& variable, double* value ) const ;
	std::string get_summary( std::string const& prefix, std::size_t target_column ) const ;

protected:
	std::string const m_missing_value ;
	genfile::SNPIdentifyingDataTest::UniquePtr m_exclusion_test ;
	typedef boost::bimap< std::string, std::pair< std::size_t, std::string > > SourceColumnMap ;
	typedef std::set< std::pair< std::string, bool > > DesiredColumns ;
	
	boost::optional< SourceColumnMap > m_column_map ;
	
	std::vector< genfile::SNPIdentifyingData2 > m_snps ;
	int m_degrees_of_freedom ;
	Eigen::MatrixXf m_betas ;
	Eigen::MatrixXf m_ses ;
	Eigen::MatrixXf m_covariance ;
	Eigen::VectorXf m_pvalues ;
	Eigen::VectorXf m_info ;
	Eigen::VectorXf m_maf ;
	Eigen::MatrixXf m_sample_counts ;
	typedef std::map< std::string, std::vector< double > > ExtraVariables ;
	ExtraVariables m_extra_variables ;
	
	FlatFileFrequentistGenomeWideAssociationResults() ;
	
	virtual void setup_columns( std::vector< std::string > const& column_names ) = 0 ;
	virtual DesiredColumns get_desired_columns() const = 0 ;
	virtual bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const = 0 ;
	virtual bool check_if_snp_accepted( std::size_t snp_index ) const = 0 ;
	virtual void store_value( int snp_index, std::string const& variable, double value ) = 0 ;
	
private:

	void setup(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) ;

	void setup(
		statfile::BuiltInTypeStatSource::UniquePtr source,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) ;
	
	SourceColumnMap get_source_column_map(
		statfile::BuiltInTypeStatSource const& source
	) const ;
	void check_columns(
		statfile::BuiltInTypeStatSource const& source,
		SourceColumnMap const& column_map
	) const ;

	void resize_storage( Eigen::MatrixXf::Index const N_snps, Eigen::MatrixXf::Index const degrees_of_freedom ) ;
	
	void free_unused_memory() ;
} ;

#endif
