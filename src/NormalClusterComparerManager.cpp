
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set.hpp"
#include "NormalClusterComparer.hpp"
#include "NormalClusterComparerManager.hpp"

void NormalClusterFitComparerManager::add_comparer( std::string const& name, NormalClusterFitComparer::UniquePtr comparer ) {
	m_comparers.insert( name, comparer ) ;
}

void NormalClusterFitComparerManager::begin_processed_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
}
void NormalClusterFitComparerManager::processed_snp( genfile::SNPIdentifyingData snp, genfile::VariantDataReader& data_reader ) {
	Eigen::MatrixXd intensities ;
	genfile::vcf::MatrixSetter< Eigen::MatrixXd > intensity_setter( intensities ) ;
	data_reader.get( m_spec, intensity_setter ) ;
	assert( std::size_t( intensities.rows() ) == 2 ) ;
	assert( std::size_t( intensities.cols() ) == m_number_of_samples ) ;
	process_snp( snp, intensities ) ;
}

void NormalClusterFitComparerManager::send_results( genfile::SNPIdentifyingData const& snp, std::string const& s1, std::string const& s2, std::string const& s3, Eigen::MatrixXd const& m ) const {
	m_results_signal( snp, s1, s2, s3, m ) ;
}

void NormalClusterFitComparerManager::process_snp( genfile::SNPIdentifyingData snp, Eigen::MatrixXd const& intensities ) const {
	for( Fits::const_iterator fit1 = m_fits.begin(); fit1 != m_fits.end(); ++fit1 ) {
		Fits::const_iterator fit2 = fit1 ;
		for( ++fit2; fit2 != m_fits.end(); ++fit2 ) {
			process_snp_fits( snp, intensities, fit1->first, fit1->second, fit2->first, fit2->second, m_comparers ) ;
		}
	}
}

void NormalClusterFitComparerManager::process_snp_fits(
	genfile::SNPIdentifyingData const& snp,
	Eigen::MatrixXd const& intensities,
	std::string const& fit1_name,
	Eigen::MatrixXd const& fit1,
	std::string const& fit2_name,
	Eigen::MatrixXd const& fit2,
	Comparers const& comparers
) const {
	for( Comparers::const_iterator i = comparers.begin(); i != comparers.end(); ++i ) {
		send_results(
			snp,
			fit1_name,
			fit2_name,
			i->first,
			i->second->compare( intensities, fit1, fit2 )
		) ;
	}
}

void NormalClusterFitComparerManager::end_processed_snps() {
	// do nothing.
}
