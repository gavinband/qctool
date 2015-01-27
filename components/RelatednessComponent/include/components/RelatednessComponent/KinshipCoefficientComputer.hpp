
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef COMPONENTS_RELATEDNESS_COMPONENT_KINSHIPCOEFFICIENTCOMPUTATION_HPP
#define COMPONENTS_RELATEDNESS_COMPONENT_KINSHIPCOEFFICIENTCOMPUTATION_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/signals2.hpp>
#include "Eigen/Core"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "worker/Worker.hpp"
#include "worker/Task.hpp"
#include "appcontext/UIContext.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/FileUtil.hpp"
#include "components/RelatednessComponent/KinshipCoefficientManager.hpp"
//#include "components/RelatednessComponent/KinshipCoefficientBlockTask.hpp"

struct KinshipCoefficientComputer: public KinshipCoefficientManager, public genfile::SNPDataSourceProcessor::Callback
{
public:
	typedef std::auto_ptr< KinshipCoefficientComputer > UniquePtr ;
	
	struct Computation: public genfile::SNPDataSourceProcessor::Callback {
		typedef Eigen::MatrixXd Matrix ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > IntegerMatrix ;
		typedef Eigen::Matrix< int, Eigen::Dynamic, 1 > IntegerVector ;
		//typedef Eigen::MatrixXf IntegerMatrix ;
		//typedef Eigen::VectorXf IntegerVector ;
		typedef std::auto_ptr< Computation > UniquePtr ;
		virtual std::size_t number_of_snps_included() const = 0 ;
		virtual Matrix const& result() const = 0 ;
		virtual IntegerMatrix const& nonmissingness() const = 0 ;
	} ;
	static genfile::SNPDataSourceProcessor::Callback::UniquePtr create_computation( std::string const& spec ) ;

public:
	~KinshipCoefficientComputer() throw() {}

	KinshipCoefficientComputer(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		appcontext::UIContext& ui_context,
		Computation::UniquePtr computation
	) ;

	void begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) ;
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader::SharedPtr data_reader ) ;
	void end_processing_snps() ;

private:
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
	std::size_t m_number_of_snps_processed ;
	genfile::CohortIndividualSource const& m_samples ;
	std::vector< Computation::Matrix > m_result ;
	std::vector< Computation::IntegerMatrix > m_non_missing_count ;
	std::vector< Computation::Vector > m_genotypes ;
	std::vector< Computation::IntegerVector > m_non_missingness ;
	Computation::UniquePtr m_computation ;
	void (*m_accumulate_xxt)( Computation::Vector*, Computation::Matrix*, int const begin_sample_i, int const end_sample_i, int const begin_sample_j, int const end_sample_j, double const ) ;
	void (*m_accumulate_xxt_integer)( Computation::IntegerVector*, Computation::IntegerMatrix*, int const begin_sample_i, int const end_sample_i, int const begin_sample_j, int const end_sample_j, double const ) ;
} ;

namespace impl {
	struct Dispatcher ;

	struct NormaliseGenotypesAndComputeXXt: public KinshipCoefficientComputer::Computation {
		static UniquePtr create( worker::Worker*, std::string const& method ) ;
		NormaliseGenotypesAndComputeXXt( worker::Worker*, std::string const& method ) ;

		Computation::Matrix const& result() const ;
		Computation::IntegerMatrix const& nonmissingness() const ;
		void begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) ;
		void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader::SharedPtr data_reader ) ;
		void end_processing_snps() ;
		std::size_t number_of_snps_included() const ;
	private:
		double const m_call_threshhold ;
		double const m_allele_frequency_threshhold ;
		std::size_t m_number_of_snps_included ;
		Computation::Matrix m_result ;
		Computation::IntegerMatrix m_nonmissingness ;
		std::auto_ptr< Dispatcher > m_dispatcher ;
	} ;
}

#endif
