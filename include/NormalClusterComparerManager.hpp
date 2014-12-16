
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_NORMAL_CLUSTER_FIT_COMPARER_MANAGER_HPP
#define QCTOOL_NORMAL_CLUSTER_FIT_COMPARER_MANAGER_HPP

#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/signals2/signal.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/VariantDataReader.hpp"
#include "ClusterFitter.hpp"
#include "NormalClusterComparer.hpp"

struct NormalClusterFitComparerManager:
	public ClusterFitter::ResultCallback,
	public genfile::SNPDataSourceProcessor::Callback
{
public:
	virtual ~NormalClusterFitComparerManager() {} ;
	NormalClusterFitComparerManager( std::string const& spec ): m_spec( spec ) {}
	void add_comparer( std::string const& name, NormalClusterFitComparer::UniquePtr comparer ) ;
	void begin_processed_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) ;
	void end_processed_snps() ;

	void send_results( genfile::SNPIdentifyingData const&, std::string const&, std::string const&, std::string const&, Eigen::MatrixXd const& ) const ;

private:
	std::string const m_spec ;
	typedef boost::ptr_map< std::string, NormalClusterFitComparer > Comparers ;
	Comparers m_comparers ;
	typedef std::map< std::string, Eigen::MatrixXd > Fits ;
	typedef Fits::const_iterator ConstFitIterator ;
	Fits m_fits ;
	
	std::size_t m_number_of_samples ;
	
	boost::signals2::signal< void ( genfile::SNPIdentifyingData const&, std::string const&, std::string const&, std::string const&, Eigen::MatrixXd const& ) > m_results_signal ;

	void process_snp( genfile::SNPIdentifyingData snp, Eigen::MatrixXd const& intensities ) const ;
	void process_snp_fits(
		genfile::SNPIdentifyingData const&,
		Eigen::MatrixXd const&,
		std::string const&,
		Eigen::MatrixXd const&,
		std::string const&,
		Eigen::MatrixXd const&,
		NormalClusterFitComparerManager::Comparers const&
	) const ;
} ;

#endif
