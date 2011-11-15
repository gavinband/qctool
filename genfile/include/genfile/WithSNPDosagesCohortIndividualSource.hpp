#ifndef GENFILE_WITHSNPDOSAGESCOHORTINDIVIDUALSOURCE_HPP
#define GENFILE_WITHSNPDOSAGESCOHORTINDIVIDUALSOURCE_HPP

#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "../config.hpp"
#if HAVE_EIGEN
#include <Eigen/Core>
#endif
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CompositeCohortIndividualSource.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	
	class SNPMatcher: public genfile::SNPIdentifyingDataTest
	{
	public:
		typedef genfile::GenomePosition GenomePosition ;
		typedef std::auto_ptr< SNPMatcher > UniquePtr ;
		typedef boost::shared_ptr< SNPMatcher > SharedPtr ;
		static UniquePtr create( std::string test_spec ) ;

		// Convenience function as we access these objects by pointer.
		bool match(
			std::string snpid,
			std::string rsid,
			GenomePosition position,
			char allele1,
			char allele2
		) const {
			return this->operator()( snpid, rsid, position, allele1, allele2 ) ;
		}
		
		virtual std::string get_spec() const = 0 ;
	public:
		std::string const m_spec ;
		virtual ~SNPMatcher() {}
		SNPMatcher() {} ;
	} ;
	
	class RSIDMatchesTest: public SNPMatcher
	{
	public:
		RSIDMatchesTest( std::string const& rsid ) ;
		bool operator()( std::string, std::string, GenomePosition position, char, char ) const ;
		std::string display() const ;
		std::string get_spec() const ;
	private:
		std::string const m_rsid ;
	} ;

	class SNPIDMatchesTest: public SNPMatcher
	{
	public:
		SNPIDMatchesTest( std::string const& snpid ) ;
		bool operator()( std::string, std::string, GenomePosition position, char, char ) const ;
		std::string display() const ;
		std::string get_spec() const ;
	private:
		std::string const m_snpid ;
	} ;

	class PositionMatchesTest: public SNPMatcher
	{
	public:
		PositionMatchesTest( GenomePosition const& position ) ;
		bool operator()( std::string, std::string, GenomePosition position, char, char ) const ;
 		std::string display() const ;
		std::string get_spec() const ;
 	private:
		GenomePosition const m_position ;
	} ;

 	class WithSNPDosagesCohortIndividualSource: public genfile::CompositeCohortIndividualSource {
	public:
		typedef std::auto_ptr< WithSNPDosagesCohortIndividualSource > UniquePtr ;
		typedef std::map< SNPMatcher::SharedPtr, std::set< std::string > > SNPDosageSpec ;
		typedef genfile::GenomePosition GenomePosition ;

		static UniquePtr create(
			CohortIndividualSource::ConstUniquePtr sample_source,
			genfile::SNPDataSource& snp_data_source,
			SNPDosageSpec const& snp_matchers
		) ;
		
	public:
		WithSNPDosagesCohortIndividualSource(
			genfile::CohortIndividualSource::ConstUniquePtr sample_source,
			genfile::SNPDataSource& snp_data_source,
			SNPDosageSpec const& snp_matchers
		) ;
		
		std::size_t get_number_of_individuals() const ;
		ColumnSpec get_column_spec() const ;

		Entry get_entry( std::size_t sample_i, std::string const& column_name ) const ;

		genfile::CohortIndividualSource const& get_parent_source() const ;
		genfile::CohortIndividualSource const& get_base_source() const ;

		std::string get_source_spec() const ;
	private:
		
		void setup( genfile::SNPDataSource& snp_data_source, SNPDosageSpec const& rsids ) ;
		void add_column( std::string const& column_name, ::Eigen::MatrixXd const& column ) ;
		
		CohortIndividualSource::ConstUniquePtr m_source ;
		std::vector< std::string > m_column_names ;
		std::vector< ColumnType > m_column_types ;
		std::vector< std::vector< Entry > > m_column_data ;
	} ;
}

#endif
