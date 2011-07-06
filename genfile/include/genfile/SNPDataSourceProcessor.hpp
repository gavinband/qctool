#ifndef GENFILE_SNP_DATA_SOURCE_PROCESSOR_HPP
#define GENFILE_SNP_DATA_SOURCE_PROCESSOR_HPP

#include <boost/function.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	class SNPDataSourceProcessor {
		// Class SNPDataSourceProcessor.
		// This class reads all of the SNPs in a SNPDataSource one by one,
		// and feeds the data out to callback objects.
	public:
		
		typedef boost::function< void ( std::size_t, std::size_t ) > ProgressCallback ;
		
		
		struct Callback {
			virtual ~Callback() ;
			virtual void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) = 0 ;
			virtual void processed_snp( SNPIdentifyingData const& , VariantDataReader& data_reader ) = 0 ;
			virtual void end_processing_snps() = 0 ;
		} ;
		
		virtual ~SNPDataSourceProcessor() ;
		virtual void add_callback( Callback& callback ) ;
		virtual void process( genfile::SNPDataSource& source, ProgressCallback = ProgressCallback() ) = 0 ;

	protected:
		virtual void call_begin_processing_snps( std::size_t const& number_of_samples, std::size_t const& number_of_snps ) const ;
		virtual void call_processed_snp(  SNPIdentifyingData const& id_data, VariantDataReader& data_reader ) const ;
		virtual void call_end_processing_snps() const ;

	private:
		std::vector< Callback* > m_callbacks ;
	} ;
	
	class SimpleSNPDataSourceProcessor: public SNPDataSourceProcessor
		// This class visits each SNP in the source one at a time.
	{
	public:
		virtual void process( genfile::SNPDataSource& source, ProgressCallback = ProgressCallback() ) ;		
	} ;
}


#endif
