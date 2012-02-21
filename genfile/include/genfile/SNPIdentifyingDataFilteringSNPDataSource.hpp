#ifndef GENFILE_SNPIDENTIFYINGDATAFILTERINGSNPDATASOURCE_HPP
#define GENFILE_SNPIDENTIFYINGDATAFILTERINGSNPDATASOURCE_HPP

#include <vector>
#include <set>
#include <memory>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/CommonSNPFilter.hpp"

namespace genfile {
	// class SNPIdentifyingDataFilteringSNPDataSource
	// This class represents a view of its parent SNPDataSource such that
	// only SNPs whose chromosome, position, id fields and alleles satisfy the
	// given test are visible.
	class SNPIdentifyingDataFilteringSNPDataSource: public SNPDataSource
	{
	public:
		typedef std::auto_ptr< SNPIdentifyingDataFilteringSNPDataSource > UniquePtr ;
		// Create a SNPIdentifyingDataFilteringSNPDataSource from the given source and the given test.
		// The test should return true for tests that are in and false for those excluded.
		static UniquePtr create( SNPDataSource::UniquePtr source, SNPIdentifyingDataTest::UniquePtr test ) ;
		
	public:
		SNPIdentifyingDataFilteringSNPDataSource( SNPDataSource::UniquePtr source, SNPIdentifyingDataTest::UniquePtr test ) ;

		operator bool() const ;
		unsigned int number_of_samples() const ;
		OptionalSnpCount total_number_of_snps() const ;
		std::string get_source_spec() const ;

		SNPDataSource const& get_parent_source() const ;
		SNPDataSource const& get_base_source() const ;
		
	private:
		void reset_to_start_impl() ;

		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;

	protected:
		SNPDataSource& source() { return *m_source ; }
		
	private:
		
		SNPDataSource::UniquePtr m_source ;
		SNPIdentifyingDataTest::UniquePtr m_test ;
	} ;
}

#endif