#ifndef QCTOOL_SAMPLE_BY_SAMPLE_COMPUTATION_HPP
#define QCTOOL_SAMPLE_BY_SAMPLE_COMPUTATION_HPP

#include <memory>
#include <vector>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"

struct SampleBySampleComputation
{
public:
	typedef std::auto_ptr< SampleBySampleComputation > UniquePtr ;
	static UniquePtr create( std::string const& name, appcontext::OptionProcessor const& options, appcontext::UIContext& ui_context ) ;
		
public:
	virtual ~SampleBySampleComputation() {}
	
	virtual void prepare(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
	) = 0 ;
	
	virtual double operator()(
		std::size_t first_sample_i,
		std::size_t second_sample_i,
		std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
	) = 0 ;
	
	virtual std::string get_summary() const = 0 ;
} ;

#endif
