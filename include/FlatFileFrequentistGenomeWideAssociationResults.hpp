
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef FLAT_FILE_FREQUENTISTGENOMEWIDEASSOCIATIONRESULTS_HPP
#define FLAT_FILE_FREQUENTISTGENOMEWIDEASSOCIATIONRESULTS_HPP

#include <boost/bimap.hpp>
#include <boost/regex.hpp>
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
		statfile::BuiltInTypeStatSource::UniquePtr source,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) ;

	std::size_t get_number_of_SNPs() const ;

	genfile::SNPIdentifyingData2 const& get_SNP( std::size_t snp_i ) const ;
	void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const ;
	void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const ;
	void get_covariance_upper_triangle( std::size_t snp_i, Eigen::VectorXd* result ) const ; 
	void get_pvalue( std::size_t snp_i, double* result ) const ;
	void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const ;
	void get_info( std::size_t snp_i, double* result ) const ;
	void get_maf( std::size_t snp_i, double* result ) const ;
	void get_frequency( std::size_t snp_i, double* result ) const ;
	void get_variable( std::size_t snp_i, std::string const& variable, std::string* value ) const ;
	std::string get_summary( std::string const& prefix, std::size_t target_column ) const ;

public:
	struct ColumnSpec {
		ColumnSpec():
			 m_required( false )
		{}
		ColumnSpec( std::string const& name, boost::regex const& pattern, std::size_t index, bool required, std::string const& type ):
			m_name( name ),
			m_simplified_name( name ),
			m_pattern( pattern ),
			m_index( index ),
			m_required( required ),
			m_type( type )
		{}
		ColumnSpec( ColumnSpec const& other ):
			m_name( other.m_name ),
			m_simplified_name( other.m_simplified_name ),
			m_pattern( other.m_pattern ),
			m_index( other.m_index ),
			m_required( other.m_required ),
			m_type( other.m_type )
		{}
		ColumnSpec& operator=( ColumnSpec const& other ) {
			m_name = other.m_name ;
			m_name = other.m_simplified_name ;
			m_pattern = other.m_pattern ;
			m_index = other.m_index ;
			m_required = other.m_required ;
			m_type = other.m_type ;
			return *this ;
		}
		bool operator<( ColumnSpec const& other ) const {
			return m_name < other.m_name ;
		}
		bool operator==( ColumnSpec const& other ) const {
			return (m_name == other.m_name) ;
		}
		std::string const& name() const { return m_name ; }
		std::string const& simplified_name() const { return m_simplified_name ; }
		boost::regex const& pattern() const { return m_pattern ; }
		std::size_t index() const { return m_index ; }
		bool const& required() const { return m_required ; }
		std::string const& type() const { return m_type ; }

		void set_simplified_name( std::string const& name ) { m_simplified_name = name ; }
	private:
		std::string m_name ;
		std::string m_simplified_name ;
		boost::regex m_pattern ;
		std::size_t m_index ;
		bool m_required ;
		std::string m_type ;
	} ;
	typedef std::vector< ColumnSpec > DesiredColumns ;
 	typedef boost::bimap< ColumnSpec, std::size_t > SourceColumnMap ;
	
protected:
	std::string const m_missing_value ;
	genfile::SNPIdentifyingDataTest::UniquePtr m_exclusion_test ;

	DesiredColumns m_desired_columns ;
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
	typedef Eigen::MatrixXf::ConstRowXpr Row ;
	typedef std::map< std::string, std::vector< std::string > > ExtraVariables ;
	ExtraVariables m_extra_variable_storage ;
	
	FlatFileFrequentistGenomeWideAssociationResults() ;
	
	virtual DesiredColumns setup_columns( std::vector< std::string > const& column_names ) = 0 ;
	virtual bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const = 0 ;
	virtual bool check_if_snp_accepted( std::size_t snp_index ) const = 0 ;
	virtual void store_value( int snp_index, std::string const& variable, std::string const& value ) = 0 ;
	
private:

	void setup(
		statfile::BuiltInTypeStatSource::UniquePtr source,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) ;
	
	SourceColumnMap get_source_column_map(
		statfile::BuiltInTypeStatSource const& source,
		DesiredColumns const&
	) const ;
	void check_columns(
		statfile::BuiltInTypeStatSource const& source,
		SourceColumnMap const& column_map
	) const ;

	void resize_storage( Eigen::MatrixXf::Index const N_snps, Eigen::MatrixXf::Index const degrees_of_freedom ) ;
	
	void free_unused_memory() ;
} ;

#endif
