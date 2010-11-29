#include "BayesFactorWeighter.hpp"
#include "BayesFactorCombiner.hpp"

struct TestWeighter: public BayesFactorWeighter::Base
{
	virtual double operator()( TranscriptDefinition const& transcript_defn, uint32_t position, double value ) const {
		if( transcript_defn.transcription_range().contains( position )) {
			return value ;
		}
		else {
			return 0
		}
	}	
} ;


struct TestCombinerContext: public BayesFactorCombinerContext
{
	TestCombinerContext( std::string const& input_string, std::vector< pathway::TranscriptDefn > const& transcript_defns )
		: m_stat_sink( m_result_stream ),
		  m_input_stream( input_string ),
  		  m_transcript_defns( transcript_defns ),
		  m_transcript_source( m_transcript_defns ),
		  m_stat_source( m_input_stream )
	) {
	}
	virtual statfile::StatSource< std::size_t, TranscriptDefinition >& transcript_source() { return m_transcript_source ; }
	virtual BayesFactorWeighter::Base const& bf_weighter() const { return m_bf_weighter ; } ;
	virtual statfile::BuiltInTypeStatSource& bf_source() const { return *m_stat_source ; }
	virtual statfile::BuiltInTypeStatSink& result_sink() const { return *m_stat_sink ; }
	virtual OstreamTee& logger() const { return m_ostream_tee ; }
	virtual uint32_t position_margin() const { return 0 ; } ;

	virtual void write_start_banner() const {} ;
	virtual void write_end_banner() const {} ;

private:
	std::ostringstream m_result_stream ;
	std::istringstream m_input_stream ;
	std::vector< pathway::TranscriptDefinition > m_transcript_defns ;
	
	TestWeighter m_bf_weighter ;
	RFormatStatSink m_stat_sink ;
	RFormatStatSource m_stat_source ;
	VectorStatSource< pathway::TranscriptDefinition > m_transcript_source ;
	OstreamTee m_ostream_tee ;
} ;
