#include <string>
#include <boost/signals2/signal.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "PairwiseCallComparer.hpp"
#include "PairwiseCallComparerManager.hpp"

namespace impl {
	struct CallComparerFileOutputter {
		typedef std::auto_ptr< CallComparerFileOutputter > UniquePtr ;
		typedef boost::shared_ptr< CallComparerFileOutputter > SharedPtr ;
		
		static UniquePtr create( std::string const& filename ) { return UniquePtr( new CallComparerFileOutputter( filename ) ) ; }

		CallComparerFileOutputter( std::string const& filename ):
			m_filename( filename ),
			m_sink( statfile::BuiltInTypeStatSink::open( filename ))
		{
			(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "callset_1" | "callset_2" | "comparison_method" | "comparison_variable" | "value" ;
		}

		void write_comparison(
			genfile::SNPIdentifyingData const& snp,
			std::string const& callset1,
			std::string const& callset2,
			std::string const& comparison_method,
			std::string const& comparison_variable,
			genfile::VariantEntry const& value
		) {
			(*m_sink)
				<< snp.get_SNPID()
				<< snp.get_rsid()
				<< std::string( snp.get_position().chromosome() )
				<< snp.get_position().position()
				<< std::string( 1, snp.get_first_allele() )
				<< std::string( 1, snp.get_second_allele() )
				<< callset1
				<< callset2
				<< comparison_method
				<< comparison_variable ;
			(*m_sink)
				<< value.as< double >()
				<< statfile::end_row() ;
			;
		}

	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	} ;
}

void PairwiseCallComparerManager::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Call comparison options" ) ;
	options[ "-compare-calls" ]
		.set_description( "Compare genotype calls from the given fields of a VCF or other file." )
		.set_takes_single_value() ;
	options[ "-compare-calls-file" ]
		.set_description( "Name of output file to put call comparisons in." )
		.set_takes_single_value()
		.set_default_value( "call-comparisons.csv" ) ;
}

PairwiseCallComparerManager::UniquePtr PairwiseCallComparerManager::create( appcontext::OptionProcessor const& options ) {
	return PairwiseCallComparerManager::create(
		options.get_value< std::string >( "-compare-calls-file" ),
		genfile::string_utils::split_and_strip_discarding_empty_entries(
			options.get_value< std::string >( "-compare-calls" ),
			",",
			" \t"
		)
	) ;
}

PairwiseCallComparerManager::UniquePtr PairwiseCallComparerManager::create( std::string const& spec, std::vector< std::string > const& call_fields ) {
	PairwiseCallComparerManager::UniquePtr result( new PairwiseCallComparerManager( call_fields ) ) ;
	result->send_results_to(
		boost::bind(
			&impl::CallComparerFileOutputter::write_comparison,
			impl::CallComparerFileOutputter::SharedPtr( impl::CallComparerFileOutputter::create( spec )),
			_1, _2, _3, _4, _5, _6
		)
	) ;
	result->add_comparer(
		"AlleleFrequencyTestCallComparer",
		PairwiseCallComparer::create( "AlleleFrequencyTestCallComparer" )
	) ;
	return result ;
}

PairwiseCallComparerManager::PairwiseCallComparerManager( std::vector< std::string > const& call_fields ):
 	m_call_fields( call_fields )
{}

void PairwiseCallComparerManager::add_comparer( std::string const& name, PairwiseCallComparer::UniquePtr comparer ) {
	m_comparers.insert( name, comparer ) ;
}

void PairwiseCallComparerManager::send_results_to( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void PairwiseCallComparerManager::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	// do nothing.
	m_calls.clear() ;
}

void PairwiseCallComparerManager::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	for( std::size_t i = 0; i < m_call_fields.size(); ++i ) {
		data_reader.get( m_call_fields[i], m_calls[ m_call_fields[i] ] ) ;
	}

	for( Calls::const_iterator i = m_calls.begin(); i != m_calls.end(); ++i ) {
		Calls::const_iterator j = i ;
		for( ++j; j != m_calls.end(); ++j ) {
			for( Comparers::const_iterator comparer_i = m_comparers.begin(); comparer_i != m_comparers.end(); ++comparer_i ) {
				std::map< std::string, genfile::VariantEntry > const& result = comparer_i->second->compare( i->second, j->second ) ;
				std::map< std::string, genfile::VariantEntry >::const_iterator
					result_i = result.begin(),
					result_end = result.end() ;
				for( ; result_i != result_end; ++result_i ) {
					m_result_signal( snp, i->first, j->first, comparer_i->first, result_i->first, result_i->second ) ;
				}
			}
		}
	}
}

void PairwiseCallComparerManager::end_processing_snps() {
	// nothing to do.
}
