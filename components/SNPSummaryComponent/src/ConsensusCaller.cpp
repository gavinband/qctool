
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/vcf/get_set.hpp"
#include "genfile/get_set_eigen.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SNPSummaryComponent/ConsensusCaller.hpp"
#include "components/SNPSummaryComponent/QuangStyleConsensusCaller.hpp"
#include "components/SNPSummaryComponent/LeastMissingConsensusCaller.hpp"

ConsensusCaller::UniquePtr ConsensusCaller::create( std::string const& model ) {
	if( model == "LeastMissing" ) {
		return ConsensusCaller::UniquePtr( new LeastMissingConsensusCaller() ) ;
	}
	else if( model == "QuangStyle" ) {
		return ConsensusCaller::UniquePtr( new QuangStyleConsensusCaller( 0.95, 1.0 ) ) ;
	} else if( model == "NoConflict" ) {
		return ConsensusCaller::UniquePtr( new QuangStyleConsensusCaller( 0.95, 0.0 ) ) ;
	} else {
		throw genfile::BadArgumentError( "ConsensusCaller::create()", "model=\"" + model + "\"" ) ;
	}
}

ConsensusCaller::SharedPtr ConsensusCaller::create_shared( std::string const& model ) {
	return ConsensusCaller::SharedPtr( ConsensusCaller::create( model ).release() ) ;
}

ConsensusCaller::ConsensusCaller()
{}

ConsensusCaller::~ConsensusCaller()
{}

void ConsensusCaller::send_results_to( ResultSignal::slot_type callback ) {
	m_result_signal.connect( callback ) ;
}

void ConsensusCaller::send_results(
	genfile::SNPIdentifyingData const& snp,
	Eigen::MatrixXd const& genotypes,
	std::map< std::string, std::vector< genfile::VariantEntry > > const& info
) {
	m_result_signal( boost::cref( snp ), boost::cref( genotypes ), boost::cref( info )) ;
}

void ConsensusCaller::begin_processing_snps( std::size_t number_of_samples ) {
	m_number_of_samples = number_of_samples ;
}

void ConsensusCaller::begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
	m_snp = snp ;
	m_call_names.clear() ;
}

void ConsensusCaller::end_comparisons() {}
