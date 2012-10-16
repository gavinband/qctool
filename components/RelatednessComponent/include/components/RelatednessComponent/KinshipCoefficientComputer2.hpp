
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef COMPONENTS_RELATEDNESS_COMPONENT_KINSHIPCOEFFICIENTCOMPUTATION2_HPP
#define COMPONENTS_RELATEDNESS_COMPONENT_KINSHIPCOEFFICIENTCOMPUTATION2_HPP

#include <vector>
#include <boost/ptr_container/ptr_deque.hpp>
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

struct KinshipCoefficientComputer2: public KinshipCoefficientManager, public genfile::SNPDataSourceProcessor::Callback
{
public:
	typedef std::auto_ptr< KinshipCoefficientComputer2 > UniquePtr ;
	struct BlockExtent ;
public:
	~KinshipCoefficientComputer2() throw() ;

	KinshipCoefficientComputer2(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker,
		appcontext::UIContext& ui_context
	) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

private:
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
	worker::Worker* m_worker ;
	Eigen::MatrixXd m_result ;
	Eigen::MatrixXd m_non_missing_count ;
	std::vector< BlockExtent > m_subdivision ;
	std::vector< Eigen::VectorXd > m_genotypes ;
	std::vector< Eigen::VectorXd > m_genotype_non_missingness ;
	boost::ptr_deque< worker::Task > m_tasks ;
	std::size_t m_number_of_snps_processed ;
	std::size_t m_data_index ;

public:
	struct BlockExtent {
		BlockExtent( int x, int y, int x_end, int y_end ):
			m_x( x ), m_y( y ), m_x_end( x_end ), m_y_end( y_end )
		{
			assert( m_x >= 0 ) ;
			assert( m_y >= 0 ) ;
			assert( m_x_end >= m_x ) ;
			assert( m_y_end >= m_y ) ;
		}
		
		BlockExtent( BlockExtent const& other ):
			m_x( other.m_x ),
			m_y( other.m_y ),
			m_x_end( other.m_x_end ),
			m_y_end( other.m_y_end )
		{}
		
		BlockExtent& operator=( BlockExtent const& other ) {
			m_x = other.m_x ;
			m_y = other.m_y ;
			m_x_end = other.m_x_end ;
			m_y_end = other.m_y_end ;
			return *this ;
		}
		
		int x() const { return m_x ; }
		int y() const { return m_y ; }
		int x_end() const { return m_x_end ; }
		int y_end() const { return m_y_end ; }
		int x_size() const { return m_x_end - m_x ; }
		int y_size() const { return m_y_end - m_y ; }
		
	private:
		int m_x, m_y, m_x_end, m_y_end ;
	} ;
} ;

#endif
