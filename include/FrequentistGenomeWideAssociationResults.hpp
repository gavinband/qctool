
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef FREQUENTISTGENOMEWIDEASSOCIATIONRESULTS_HPP
#define FREQUENTISTGENOMEWIDEASSOCIATIONRESULTS_HPP

#include <memory>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData2.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/VariantEntry.hpp"

struct FrequentistGenomeWideAssociationResults: public boost::noncopyable {
public:
	typedef std::auto_ptr< FrequentistGenomeWideAssociationResults > UniquePtr ;
	typedef boost::function< void ( std::size_t i, genfile::SNPIdentifyingData2 const& snp ) > GetSNPCallback ;
	typedef
		boost::function< void ( genfile::SNPIdentifyingData2 const&, std::string const&, genfile::VariantEntry const& ) > 
		SNPResultCallback ;
	typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;

	static UniquePtr create(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		boost::optional< std::string > const& effect_size_column_regex,
		std::vector< std::string > const& columns,
		genfile::SNPIdentifyingDataTest::UniquePtr test,
		SNPResultCallback callback = SNPResultCallback(),
		ProgressCallback progress_callback = ProgressCallback()
	) ;
	
public:
	virtual ~FrequentistGenomeWideAssociationResults() {}
	
	virtual void add_variable( std::string const& variable ) = 0 ;
	virtual void set_effect_size_column_regex( std::string const& ) = 0 ;
	virtual std::size_t get_number_of_SNPs() const = 0 ;
	virtual genfile::SNPIdentifyingData2 const& get_SNP( std::size_t snp_i ) const = 0 ;
	virtual int get_number_of_effect_parameters() const = 0 ;
	virtual void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ;
	virtual void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ; 
	virtual void get_covariance_upper_triangle( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ; 
	virtual void get_pvalue( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ;
	virtual void get_info( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_maf( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_frequency( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_variable( std::size_t snp_i, std::string const& variable, double* result ) const = 0 ;
	
	virtual std::string get_summary( std::string const& prefix = "", std::size_t target_column = 80 ) const = 0 ;

	// Legacy function: this should probably be removed.
	void get_SNPs( GetSNPCallback callback ) const {
		for( std::size_t i = 0; i < get_number_of_SNPs(); ++i ) {
			callback( i, get_SNP( i ) ) ;
		}
	}
} ;

#endif
