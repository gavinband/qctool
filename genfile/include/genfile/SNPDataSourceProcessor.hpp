
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNP_DATA_SOURCE_PROCESSOR_HPP
#define GENFILE_SNP_DATA_SOURCE_PROCESSOR_HPP

#include <boost/function.hpp>
#include <boost/optional.hpp>
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
		
		typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;
		
		struct Callback {
			typedef std::auto_ptr< Callback > UniquePtr ;
			virtual ~Callback() ;
			virtual void begin_processing_snps( std::size_t number_of_samples ) = 0 ;
			virtual void processed_snp( SNPIdentifyingData const&, VariantDataReader& data_reader ) = 0 ;
			virtual void end_processing_snps() = 0 ;
		} ;
		
		virtual ~SNPDataSourceProcessor() ;
		virtual void add_callback( Callback::UniquePtr callback ) ;
		virtual void add_callback( Callback& callback ) ;
		virtual void process( genfile::SNPDataSource& source, ProgressCallback = ProgressCallback() ) = 0 ;

		std::vector< Callback* > const& get_callbacks() const { return m_callbacks ; }

	protected:
		virtual void call_begin_processing_snps( std::size_t const& number_of_samples ) const ;
		virtual void call_processed_snp(  SNPIdentifyingData const& id_data, VariantDataReader& data_reader ) const ;
		virtual void call_end_processing_snps() const ;

	private:
		std::vector< Callback* > m_callbacks ;
		std::vector< bool > m_should_manage ;
	} ;
	
	class SimpleSNPDataSourceProcessor: public SNPDataSourceProcessor
		// This class visits each SNP in the source one at a time.
	{
	public:
		virtual void process( genfile::SNPDataSource& source, ProgressCallback = ProgressCallback() ) ;		
	} ;
}


#endif
