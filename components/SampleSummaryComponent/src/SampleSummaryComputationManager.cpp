#include <string>
#include <memory>
#include <boost/noncopyable.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputationManager.hpp"

SampleSummaryComputationManager::UniquePtr SampleSummaryComputationManager::create() {
	return SampleSummaryComputationManager::UniquePtr( new SampleSummaryComputationManager() ) ;
}

void SampleSummaryComputationManager::add( std::string const& name, std::string const& chromosome_spec, SampleSummaryComputation::UniquePtr computation ) {
	genfile::Chromosome chr( chromosome_spec ) ;
	if( chr == genfile::Chromosome() ) {
		if( chromosome_spec != "all chromosomes" && chromosome_spec != "autosomal chromosomes" && chromosome_spec != "sex chromosomes" ) {
			throw genfile::BadArgumentError( "SampleSummaryComputationManager::add()", "chromosome_spec=\"" + chromosome_spec + "\"" ) ;
		}
	}
	m_computations.insert( std::make_pair( name, chromosome_spec ), computation ) ;
}

void SampleSummaryComputationManager::begin_processing_snps( std::size_t number_of_samples ) {
	m_snp_index = 0 ;
	m_genotypes.resize( number_of_samples, 3 ) ;
}

void SampleSummaryComputationManager::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	{
		genfile::vcf::GenotypeSetter< Eigen::MatrixBase< SampleSummaryComputation::Genotypes > > setter( m_genotypes ) ;
		data_reader.get( "genotypes", setter ) ;
	}
	Computations::iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		genfile::Chromosome chr( i->first.second ) ;
		if(
			( chr != genfile::Chromosome() && chr == snp.get_position().chromosome() )
			|| ( i->first.second == "all chromosomes" )
			|| ( i->first.second == "autosomal chromosomes" && snp.get_position().chromosome().is_autosome() )
			|| ( i->first.second == "sex chromosomes" && snp.get_position().chromosome().is_sex_determining() )
		) {
			i->second->accumulate(
				snp,
				m_genotypes,
				data_reader
			) ;
		}
	}
	++m_snp_index ;
}

void SampleSummaryComputationManager::end_processing_snps() {
	Computations::iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		i->second->compute(
			boost::bind(
				boost::ref( m_result_signal ),
				i->first.first,
				_1,
				_2,
				i->first.first + " for " + i->first.second,
				_3
			)
		) ;
	}
}

void SampleSummaryComputationManager::add_result_callback(  ResultSignal::slot_type callback ) {
	m_result_signal.connect( callback ) ;
}

std::string SampleSummaryComputationManager::get_summary( std::string const& prefix, std::size_t column_width ) {
	std::string result ;
	Computations::const_iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		result += i->second->get_summary( prefix, column_width ) + "\n";
	}
	return result ;
}
