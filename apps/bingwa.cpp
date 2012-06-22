
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <sstream>
#include <iomanip>
#include <memory>
#include <set>
#include <map>
#include <utility>
#include <boost/bimap.hpp>
#include <boost/bind.hpp>
#include <boost/signals2/signal.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/regex.hpp>
#include <boost/math/distributions/normal.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/ProgramFlow.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/utility.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPIdentifyingData2.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/CommonSNPFilter.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"

namespace globals {
	std::string const program_name = "bingwa" ;
	std::string const program_version = "0.2" ;
}

struct FrequentistGenomeWideAssociationResults: public boost::noncopyable {
	typedef std::auto_ptr< FrequentistGenomeWideAssociationResults > UniquePtr ;
	virtual ~FrequentistGenomeWideAssociationResults() {}
	typedef boost::function< void ( std::size_t i, genfile::SNPIdentifyingData2 const& snp ) > SNPCallback ;
	virtual std::size_t get_number_of_SNPs() const = 0 ;
	virtual genfile::SNPIdentifyingData2 const& get_SNP( std::size_t snp_i ) const = 0 ;
	virtual void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ;
	virtual void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ; 
	virtual void get_pvalue( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ;
	virtual void get_info( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_maf_in_controls( std::size_t snp_i, double* result ) const = 0 ;
	virtual std::string get_summary( std::string const& prefix = "", std::size_t target_column = 80 ) const = 0;

	// Legacy function: this should be removed.
	void get_SNPs( SNPCallback callback ) const {
		for( std::size_t i = 0; i < get_number_of_SNPs(); ++i ) {
			callback( i, get_SNP( i ) ) ;
		}
	}
} ;

namespace {
	struct SNPExclusionTest: public boost::noncopyable {
		typedef std::auto_ptr< SNPExclusionTest > UniquePtr ;
		virtual ~SNPExclusionTest() {} ;
		virtual bool operator()( genfile::SNPIdentifyingData const& snp, double const maf, double const info ) const = 0 ;
		virtual std::string get_summary( std::string const& prefix, std::size_t target_column ) const = 0 ;
	} ;

	struct SNPIdentifyingDataExclusionTest: public SNPExclusionTest {
		SNPIdentifyingDataExclusionTest( genfile::SNPIdentifyingDataTest::UniquePtr test ):
			m_test( test )
		{
			assert( m_test.get() ) ;
		}

		bool operator()( genfile::SNPIdentifyingData const& snp, double const maf, double const info ) const {
			return (*m_test)( snp ) ;
		}
		
		std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
			return prefix + m_test->display() ;
		}
		
	private:
		genfile::SNPIdentifyingDataTest::UniquePtr m_test ;
	} ;

	struct SNPExclusionTestConjunction: public SNPExclusionTest {
		typedef std::auto_ptr< SNPExclusionTestConjunction > UniquePtr ;
		
		void add_subtest( SNPExclusionTest::UniquePtr subtest ) {
			assert( subtest.get() ) ;
			m_subtests.push_back( subtest ) ;
		}

		void add_subtest( genfile::SNPIdentifyingDataTest::UniquePtr test ) {
			assert( test.get() ) ;
			SNPExclusionTest::UniquePtr subtest( new SNPIdentifyingDataExclusionTest( test ) ) ;
			m_subtests.push_back( subtest ) ;
		}

		bool operator()( genfile::SNPIdentifyingData const& snp, double const maf, double const info ) const {
			for( std::size_t i = 0; i < m_subtests.size(); ++i ) {
				if( !m_subtests[i]( snp, maf, info )) {
					return false ;
				}
			}
			return true ;
		}
		
		std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
			std::string result = prefix + "( ";
			for( std::size_t i = 0; i < m_subtests.size(); ++i ) {
				if( i > 0 ) {
					result += " ) AND ( " ;
				}
				result += m_subtests[i].get_summary( "", target_column ) ;
			}
			result += " )" ;
			return result ;
		}
		
	private:
		boost::ptr_vector< SNPExclusionTest > m_subtests ;
	} ;

	struct ThreshholdingSummaryStatisticTest: public SNPExclusionTest {
	public:
		typedef std::auto_ptr< ThreshholdingSummaryStatisticTest > UniquePtr ;
	private:
		typedef std::map< std::string, std::pair< double, double > > VariableBounds ;
	public:
		void set_inclusion_bounds( std::string const& variable, double lower_bound, double upper_bound ) {
			if( variable != "controls_maf" && variable != "info" ) {
				throw genfile::BadArgumentError( "ThreshholdingSummaryStatisticTest::set_inclusion_bounds()", "variable=\"" + variable + "\"" ) ;
			}
			m_bounds[ variable ] = std::make_pair( lower_bound, upper_bound ) ;
		}
		
		bool operator()( genfile::SNPIdentifyingData const& snp, double const maf, double const info ) const {
			VariableBounds::const_iterator i = m_bounds.begin(), end_i = m_bounds.end() ;
			double value ;
			for( ; i != end_i; ++i ) {
				if( i->first == "controls_maf" ) {
					value = maf ;
				}
				else if( i->first == "info" ) {
					value = info ;
				}
				else {
					assert(0) ;
				}
				if( value < i->second.first || value > i->second.second ) {
					return false ;
				}
			}
			return true ;
		}

		std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
			std::string result = prefix + "( " ;
			VariableBounds::const_iterator i = m_bounds.begin(), end_i = m_bounds.end() ;
			using genfile::string_utils::to_string ;
			for( std::size_t count = 0 ; i != end_i; ++i, ++count ) {
				if( count > 0 ) {
					result += " ) AND ( " ;
				}
				result += i->first + " in [" + to_string( i->second.first ) + ", " + to_string( i->second.second ) + "]" ;
			}
			result += " )" ;
			return result ;
		}

	private:
		VariableBounds m_bounds ;
	} ;
}

struct SNPTESTResults: public FrequentistGenomeWideAssociationResults {
	typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;
	typedef
		boost::function< void ( std::size_t, genfile::SNPIdentifyingData2 const&, std::string const&, std::string const&, genfile::VariantEntry const& ) > 
		SNPResultCallback ;
	
	SNPTESTResults(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		SNPExclusionTest::UniquePtr test,
		SNPResultCallback callback = SNPResultCallback(),
		ProgressCallback progress_callback = ProgressCallback()
	):
		m_filenames( filenames ),
	 	m_exclusion_test( test )
	{
		setup( filenames, callback, progress_callback ) ;
	}

	std::size_t get_number_of_SNPs() const {
		return m_snps.size() ;
	}

	genfile::SNPIdentifyingData2 const& get_SNP( std::size_t snp_i ) const {
		assert( snp_i < m_snps.size() ) ;
		return m_snps[ snp_i ] ;
	}

	void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const {
		*result = m_betas.row( snp_i ).cast< double >() ;
	}
	void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const {
		*result = m_ses.row( snp_i ).cast< double >() ;
	}
	void get_pvalue( std::size_t snp_i, double* result ) const {
		*result = m_pvalues( snp_i ) ;
	}
	void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const {
		*result = m_sample_counts.row( snp_i ).cast< double >() ;
	}
	void get_info( std::size_t snp_i, double* result ) const {
		*result = m_info( snp_i ) ;
	}
	void get_maf_in_controls( std::size_t snp_i, double* result ) const {
		*result = m_controls_maf( snp_i ) ;
	}

	std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
		using genfile::string_utils::to_string ;
		
		// estimate memory used in SNPs.
		unsigned long mem_used = 0 ;
		for( std::size_t i = 0; i < m_snps.size(); ++i ) {
			mem_used += m_snps[i].get_estimated_bytes_used() ;
		}
		
		std::string result = "SNPTESTResults object ("
			+ to_string( m_snps.size() )
			+ " SNPs"
			+ ", ~"
			+ to_string(
				(
					(
						m_betas.size()
						+ m_ses.size()
						+ m_pvalues.size()
						+ m_info.size()
						+ m_controls_maf.size()
						+ m_sample_counts.size()
					) * sizeof( float )
					+ mem_used
				)
				/ 1000000.0
			) + "Mb in use.)" ;
		return result ;
	}

private:
	std::vector< genfile::wildcard::FilenameMatch > const m_filenames ;
	SNPExclusionTest::UniquePtr m_exclusion_test ;
	std::vector< genfile::SNPIdentifyingData2 > m_snps ;
	Eigen::MatrixXf m_betas ;
	Eigen::MatrixXf m_ses ;
	Eigen::VectorXf m_pvalues ;
	Eigen::VectorXf m_info ;
	Eigen::VectorXf m_controls_maf ;
	Eigen::MatrixXf m_sample_counts ;
	typedef boost::bimap< std::string, std::size_t > ColumnMap ;
	
	void setup(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) {
		progress_callback( 0, filenames.size() * 100 ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			setup( i, filenames.size(), filenames[i], callback, progress_callback ) ;
		}
	}
	
	void setup(
		std::size_t const filename_i,
		std::size_t const number_of_files,
		genfile::wildcard::FilenameMatch const& filename,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) {
		statfile::BuiltInTypeStatSource::UniquePtr source( statfile::BuiltInTypeStatSource::open( filename.filename() )) ;

		ColumnMap column_map = get_columns_to_store( *source ) ;
		
		int degrees_of_freedom = ( column_map.left.find( "_beta_2" ) == column_map.left.end() ) ? 1 : 2 ;
		
		m_snps.resize( source->number_of_rows() ) ;
		m_betas.resize( source->number_of_rows(), degrees_of_freedom ) ;
		m_ses.resize( source->number_of_rows(), degrees_of_freedom ) ;
		m_pvalues.resize( source->number_of_rows() ) ;
		m_info.resize( source->number_of_rows() ) ;
		m_controls_maf.resize( source->number_of_rows() ) ;
		m_sample_counts.resize( source->number_of_rows(), 4 ) ;

		genfile::SNPIdentifyingData snp ;
		Eigen::VectorXd betas ;
		Eigen::VectorXd ses ;
		
		std::size_t snp_index = 0 ;
		for(
			;
			(*source) >> snp.SNPID() >> snp.rsid() >> snp.position().chromosome() >> snp.position().position() >> snp.first_allele() >> snp.second_allele();
			(*source) >> statfile::ignore_all()
		) {
			// Just read the values we need
			ColumnMap::right_const_iterator i = column_map.right.begin(), end_i = column_map.right.end() ;
			double value ;
			
			for( ; i != end_i; ++i ) {
				(*source)
					>> statfile::ignore( i->first - source->current_column() )
					>> value ;
				store_value( snp_index, i->second, value ) ;
			}

			if( !m_exclusion_test.get() || m_exclusion_test->operator()( snp, m_controls_maf( snp_index ), m_info( snp_index ))) {
				m_snps[ snp_index++ ] = snp ;
			}

			if( progress_callback ) {
				progress_callback( 100 * ( filename_i + double( source->number_of_rows_read() + 1 ) / source->number_of_rows() ), 100 * number_of_files ) ;
			}
		}
		
		// Now deallocate any unused memory.
		{
			m_snps.resize( snp_index ) ;
			{
				std::vector< genfile::SNPIdentifyingData2 > snps( m_snps.begin(), m_snps.end() ) ;
				m_snps.swap( snps ) ;
			}
			{
				Eigen::MatrixXf betas = m_betas.block( 0, 0, snp_index, m_betas.cols() ) ;
				m_betas.swap( betas ) ;
			}
			{
				Eigen::MatrixXf ses = m_ses.block( 0, 0, snp_index, m_betas.cols() ) ; ;
				m_ses.swap( ses ) ;
			}
			{
				Eigen::VectorXf pvalues = m_pvalues.head( snp_index ) ;
				m_pvalues.swap( pvalues ) ;
			}
			{
				Eigen::VectorXf info = m_info.head( snp_index ) ;
				m_info.swap( info ) ;
			}
			{
				Eigen::VectorXf controls_maf = m_controls_maf.head( snp_index ) ;
				m_controls_maf.swap( controls_maf ) ;
			}
			{
				Eigen::MatrixXf sample_counts = m_sample_counts.block( 0, 0, snp_index, m_sample_counts.cols() ) ;
				m_sample_counts.swap( sample_counts ) ;
			}
		}
	}
	
	ColumnMap get_columns_to_store(
		statfile::BuiltInTypeStatSource const& source
	) {
		std::set< std::string > desired_columns ;
		desired_columns.insert( "_beta_1" ) ;
		desired_columns.insert( "_beta_2" ) ;
		desired_columns.insert( "_se_1" ) ;
		desired_columns.insert( "_se_2" ) ;
		desired_columns.insert( "_pvalue" ) ;
		desired_columns.insert( "info" ) ;
		desired_columns.insert( "controls_maf" ) ;
		desired_columns.insert( "all_AA" ) ;
		desired_columns.insert( "all_AB" ) ;
		desired_columns.insert( "all_BB" ) ;
		desired_columns.insert( "all_NULL" ) ;

		ColumnMap result ;
		for( std::size_t i = 0; i < source.number_of_columns(); ++i ) {
			std::string name = source.name_of_column( i ) ;
			for( std::set< std::string >::iterator j = desired_columns.begin(); j != desired_columns.end(); ++j ) {
				if( name.size() >= j->size() && name.compare( name.size() - j->size(), j->size(), *j ) == 0 ) {
					result.insert( ColumnMap::value_type( *j, i )) ;
					//std::cerr << "Column \"" << j->first << "\" found at position " << i << ".\n" ;
				}
			}
		}
		
		std::set< std::string > required_columns = desired_columns ;
		required_columns.erase( "_beta_2" ) ;
		required_columns.erase( "_se_2" ) ;
		
		for( std::set< std::string >::const_iterator i = required_columns.begin(); i != required_columns.end(); ++i ) {
			if( result.left.find( *i ) == result.left.end() ) {
				throw genfile::MalformedInputError( source.get_source_spec(), 0 ) ;
			}
		}

		return result ;
	}
	
	void store_value(
		int snp_index,
		std::string const& variable,
		double const value
	) {
		using genfile::string_utils::to_repr ;
		
		if( variable == "_beta_1" ) {
			m_betas( snp_index, 0 ) = value ;
		}
		else if( variable == "_se_1" ) {
			m_ses( snp_index, 0 ) = value ;
		}
		if( variable == "_beta_2" ) {
			m_betas( snp_index, 1 ) = value ;
		}
		else if( variable == "_se_2" ) {
			m_ses( snp_index, 1 ) = value ;
		}
		else if( variable == "_pvalue" ) {
			m_pvalues( snp_index ) = value ;
		}
		else if( variable == "info" ) {
			m_info( snp_index ) = value ;
		}
		else if( variable == "controls_maf" ) {
			m_controls_maf( snp_index ) = value ;
		}
		else if( variable == "all_AA" ) {
			m_sample_counts( snp_index, 0 ) = value ;
		}
		else if( variable == "all_AB" ) {
			m_sample_counts( snp_index, 1 ) = value ;
		}
		else if( variable == "all_BB" ) {
			m_sample_counts( snp_index, 2 ) = value ;
		}
		else if( variable == "all_NULL" ) {
			m_sample_counts( snp_index, 3 ) = value ;
		}
	}
} ;

#if 0
struct FilteringGenomeWideAssociationResults: public FrequentistGenomeWideAssociationResults {
	FilteringGenomeWideAssociationResults(
		FrequentistGenomeWideAssociationResults::UniquePtr source,
		SNPExclusionTest::UniquePtr test
	):
		m_source( source ),
		m_test( test )
	{
		assert( m_test.get() ) ;
		link_indices() ;
	}

	// typedef boost::function< void ( std::size_t i, genfile::SNPIdentifyingData2 const& snp ) > SNPCallback ;

	std::size_t get_number_of_SNPs() const {
		return m_links.size() ;
	}
	
	genfile::SNPIdentifyingData2 const& get_SNP( std::size_t snp_i ) const {
		assert( snp_i < m_links.size() ) ;
		return m_source->get_SNP( m_links[ snp_i ] ) ;
	}

	void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const {
		assert( snp_i < m_links.size() ) ;
		return m_source->get_betas( m_links[ snp_i ], result ) ;
	} ;
	void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const {
		assert( snp_i < m_links.size() ) ;
		return m_source->get_ses( m_links[ snp_i ], result ) ;
	}
	void get_pvalue( std::size_t snp_i, double* result ) const {
		assert( snp_i < m_links.size() ) ;
		return m_source->get_pvalue( m_links[ snp_i ], result ) ;
	}
	void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const {
		assert( snp_i < m_links.size() ) ;
		return m_source->get_counts( m_links[ snp_i ], result ) ;
	}
	void get_info( std::size_t snp_i, double* result ) const {
		assert( snp_i < m_links.size() ) ;
		return m_source->get_info( m_links[ snp_i ], result ) ;
	}
	void get_maf_in_controls( std::size_t snp_i, double* result ) const {
		assert( snp_i < m_links.size() ) ;
		return m_source->get_maf_in_controls( m_links[ snp_i ], result ) ;
	}
	std::string get_summary( std::string const& prefix = "", std::size_t target_column = 80 ) const {
		std::string result
			= prefix
			+ "FilteringGenomeWideAssociationResults ("
			+ genfile::string_utils::to_string( m_links.size() )
			+ " SNPs, filtered from "
			+ m_source->get_summary( "", target_column )
			+ ".\n" ;
		result
			+= prefix
			+ "   using SNP filter "
			+ m_test->get_summary( "", 100 )
			+ "." ;
		return result ;
	}
	
private:
	FrequentistGenomeWideAssociationResults::UniquePtr m_source ;
	SNPExclusionTest::UniquePtr m_test ;
	std::vector< std::size_t > m_links ;
private:
	
	// populate the link between our filtered SNP indices and the indices in the base source.
	void link_indices() {
		std::size_t const N = m_source->get_number_of_SNPs() ;
		m_links.reserve( N ) ;
		std::size_t m_i = 0 ;
		for( std::size_t snp_i = 0; snp_i < N; ++snp_i ) {
			if( (*m_test)( *m_source, snp_i ) ) {
				m_links.resize( m_i + 1 ) ;
				m_links[ m_i ] = snp_i ;
				++m_i ;
			}
		}
	}
} ;
#endif

struct AmetComputation: public boost::noncopyable {
	typedef std::auto_ptr< AmetComputation > UniquePtr ;
	static UniquePtr create( std::string const& name, appcontext::OptionProcessor const& ) ;
	virtual ~AmetComputation() {}
	typedef genfile::SNPIdentifyingData2 SNPIdentifyingData2 ;
	typedef boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;

	struct DataGetter: public boost::noncopyable {
		virtual ~DataGetter() {} ;
		virtual std::size_t get_number_of_cohorts() const = 0 ;
		virtual bool is_non_missing( std::size_t i ) const = 0 ;
		virtual void get_counts( std::size_t, Eigen::VectorXd* result ) const = 0 ;
		virtual void get_betas( std::size_t i, Eigen::VectorXd* result ) const = 0 ;
		virtual void get_ses( std::size_t i, Eigen::VectorXd* result  ) const = 0 ;
		virtual void get_pvalue( std::size_t i, double* result ) const = 0 ;
		virtual void get_info( std::size_t i, double* result ) const = 0 ;
	} ;

	virtual void operator()(
		SNPIdentifyingData2 const&,
		DataGetter const& data_getter,
		ResultCallback callback
	) = 0 ;
	virtual std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const = 0 ;
	virtual std::string get_spec() const = 0 ;
} ;

struct PerCohortValueReporter: public AmetComputation {
	void operator()(
		SNPIdentifyingData2 const&,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = data_getter.get_number_of_cohorts() ;
		for( std::size_t i = 0; i < N; ++i ) {
			if( data_getter.is_non_missing( i ) ) {
				Eigen::VectorXd betas ;
				Eigen::VectorXd ses ;
				Eigen::VectorXd counts ;
				double pvalue ;
				double info ;
				data_getter.get_betas( i, &betas ) ;
				data_getter.get_ses( i, &ses ) ;
				data_getter.get_counts( i, &counts ) ;
				data_getter.get_pvalue( i, &pvalue ) ;
				data_getter.get_info( i, &info ) ;

				assert( counts.size() == 4 ) ;
				using genfile::string_utils::to_string ;
				std::string prefix = "cohort_" + to_string( i + 1 ) + "/";
				callback( prefix + "AA", counts(0) ) ;
				callback( prefix + "AB", counts(1) ) ;
				callback( prefix + "BB", counts(2) ) ;
				callback( prefix + "NULL", counts(3) ) ;
				
				assert( betas.size() == ses.size() ) ;
				for( int j = 0; j < betas.size(); ++j ) {
					callback( prefix + "beta_" + to_string( j+1 ), betas(j) ) ;
					callback( prefix + "se_" + to_string( j+1 ), betas(j) ) ;
				}
				callback( prefix + "pvalue", pvalue ) ;
				callback( prefix + "info", info ) ;
			}
		}
	}
	
	std::string get_spec() const {
		return "PerCohortValueReporter" ;
	}
	
	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}
} ;

namespace impl {
	bool get_betas_and_ses_one_per_study( AmetComputation::DataGetter const& data_getter, Eigen::VectorXd& betas, Eigen::VectorXd& ses, Eigen::VectorXd& non_missingness ) {
		std::size_t N = betas.size() ;
		for( std::size_t i = 0; i < N; ++i ) {
			non_missingness( i ) = data_getter.is_non_missing( i ) ? 1.0 : 0.0 ;
			if( non_missingness( i )) {
				Eigen::VectorXd data ;
				data_getter.get_betas( i, &data ) ;
				if( data.size() != 1 ) {
					return false ;
				}
				betas(i) = data(0) ;
				data_getter.get_ses( i, &data ) ;
				if( data.size() != 1 ) {
					return false ;
				}
				ses(i) = data(0) ;
				
				// check that betas and ses are 

				if( betas(i) != betas(i) || ses(i) != ses(i) ) {
					non_missingness(i) = 0 ;
				}
			}
		}
		return true ;
	}
	
	std::string serialise( Eigen::VectorXd const& vector ) {
		std::ostringstream ostr ;
		for( int i = 0; i < vector.size(); ++i ) {
			if( i > 0 ) {
				ostr << "," ;
			}
			if( vector(i) == vector(i) ) {
				ostr << vector(i) ;
			} else {
				ostr << "NA" ;
			}
		}
		return ostr.str() ;
	}
}

struct FixedEffectFrequentistMetaAnalysis: public AmetComputation {
	void operator()(
		SNPIdentifyingData2 const& snp,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = data_getter.get_number_of_cohorts() ;
		if( N == 0 ) {
			return ;
		}

		Eigen::VectorXd betas( N ) ;
		Eigen::VectorXd ses( N ) ;
		Eigen::VectorXd non_missingness( N ) ;
		
		if( !impl::get_betas_and_ses_one_per_study( data_getter, betas, ses, non_missingness ) ) {
			callback( "FixedEffectFrequentistMetaAnalysis/meta_beta", genfile::MissingValue() ) ;
			callback( "FixedEffectFrequentistMetaAnalysis/meta_se", genfile::MissingValue() ) ;
			callback( "FixedEffectFrequentistMetaAnalysis/meta_pvalue", genfile::MissingValue() ) ;
			return ;
		}
		else {
			Eigen::VectorXd inverse_variances = ( ses.array() * ses.array() ).inverse() ;
			for( int i = 0; i < int(N); ++i ) {
				if( non_missingness( i ) == 0.0 ) {
					inverse_variances( i ) = 0 ;
					betas( i ) = 0 ;
					ses( i ) = 0 ;
				}
			}
			double const meta_beta = ( inverse_variances.array() * betas.array() ).sum() / inverse_variances.sum() ;
			double const meta_se = std::sqrt( 1.0 / inverse_variances.sum() ) ;

			callback( "FixedEffectFrequentistMetaAnalysis/meta_beta", meta_beta ) ;
			callback( "FixedEffectFrequentistMetaAnalysis/meta_se", meta_se ) ;

			//std::cerr << "SNP: " << snp << ": betas = " << betas << ", ses = " << ses << ".\n" ;

			if( meta_se > 0 && meta_se != std::numeric_limits< double >::infinity() ) {
				typedef boost::math::normal NormalDistribution ;
				NormalDistribution normal( 0, meta_se ) ;
				// P-value is the mass under both tails of the normal distribution larger than |meta_beta|
				double const pvalue = 2.0 * boost::math::cdf( boost::math::complement( normal, std::abs( meta_beta ) ) ) ;
				callback( "FixedEffectFrequentistMetaAnalysis/pvalue", pvalue ) ;
			}
			
		}
	}

	std::string get_spec() const {
		return "FixedEffectFrequentistMetaAnalysis" ;
	}

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}
} ;

//
// Bayesian meta-analysis treating the prior as normal with mean 0 and variance matrix \Sigma
// and the likelihood as normal and given by
// the estimated betas and ses.
struct ApproximateBayesianMetaAnalysis: public AmetComputation {
	ApproximateBayesianMetaAnalysis(
		std::string const& name,
		Eigen::MatrixXd const& sigma
	):
		m_name( name ),
		m_prefix( "ApproximateBayesianMetaAnalysis/" + name ),
		m_sigma( sigma )
	{
		assert( m_sigma.rows() == m_sigma.cols() ) ;
	}

	void operator()(
		SNPIdentifyingData2 const& snp,
		DataGetter const& data_getter,
		ResultCallback callback
	) {
		std::size_t const N = data_getter.get_number_of_cohorts() ;
		if( N == 0 ) {
			return ;
		}

		Eigen::VectorXd betas( N ) ;
		Eigen::VectorXd ses( N ) ;
		Eigen::VectorXd non_missingness( N ) ;
		
		if(
			!impl::get_betas_and_ses_one_per_study( data_getter, betas, ses, non_missingness )
			|| non_missingness.sum() == 0
		) {
			callback( m_prefix + "bf", genfile::MissingValue() ) ;
			return ;
		}
		else {
			// deal with missingness.
			// Make a matrix that will select the rows and columns we want.
			assert( betas.size() == m_sigma.rows() ) ;
			
			Eigen::MatrixXd prior_selector = Eigen::MatrixXd::Zero( non_missingness.sum(), betas.size() ) ;
			Eigen::MatrixXd prior ;
			{
				int count = 0 ;
				for( int i = 0; i < betas.size(); ++i ) {
					if( non_missingness(i) ) {
						prior_selector( count++, i ) = 1 ;
					}
				}
				
				prior = prior_selector * m_sigma * prior_selector.transpose() ;
				betas = prior_selector * betas ;
				ses = prior_selector * ses ;
			}
			
			assert( ses.size() == non_missingness.sum() ) ;
			assert( betas.size() == non_missingness.sum() ) ;
			assert( prior.rows() == non_missingness.sum() ) ;
			assert( prior.cols() == non_missingness.sum() ) ;
			
			Eigen::MatrixXd V = ses.asDiagonal() ;
			V = V.array().square() ;
			
			// std::cerr << std::resetiosflags( std::ios::floatfield ) ;
			// std::cerr << "V = " << V << ".\n" ;
			// std::cerr << "prior = " << prior << ".\n" ;
			// std::cerr << "betas = " << betas.transpose() << ".\n" ;
			// std::cerr << "ses = " << ses.transpose() << ".\n" ;
			
			// I hope LDLT copes with noninvertible matrices.
			// Maybe it doesn't...but let's find out.
			Eigen::LDLT< Eigen::MatrixXd > Vsolver( V ) ;
			Eigen::LDLT< Eigen::MatrixXd > V_plus_prior_solver( V + prior ) ;
			
			Eigen::VectorXd exponent = betas.transpose() * ( Vsolver.solve( betas ) - V_plus_prior_solver.solve( betas ) ) ;
			
			assert( exponent.size() == 1 ) ;
			
			double const constant = std::sqrt( Vsolver.vectorD().prod() / V_plus_prior_solver.vectorD().prod() ) ;
			double const result = constant * std::exp( 0.5 * exponent(0) ) ;
			
			callback( m_prefix + "/bf", result ) ;
		}
	}

	std::string get_spec() const {
		return "ApproximateBayesianMetaAnalysis( " + m_name + " ) with prior:\n" + genfile::string_utils::to_string( m_sigma ) ;
	}

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}
private:
	std::string const m_name ;
	std::string const m_prefix ;
	Eigen::MatrixXd const m_sigma ;
} ;


AmetComputation::UniquePtr AmetComputation::create( std::string const& name, appcontext::OptionProcessor const& options ) {
	AmetComputation::UniquePtr result ;
	if( name == "FixedEffectFrequentistMetaAnalysis" ) {
		result.reset( new FixedEffectFrequentistMetaAnalysis() ) ;
	}
	else if( name == "PerCohortValueReporter" ) {
		result.reset( new PerCohortValueReporter() ) ;
	}
	else {
		throw genfile::BadArgumentError( "AmetComputation::create()", "name=\"" + name + "\"" ) ;
	}
	return result ;
}

struct AmetOptions: public appcontext::CmdLineOptionProcessor {
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;

		options.declare_group( "File handling options" ) ;
		options[ "-snptest" ]
			.set_description( "Specify the path of a file containing SNPTEST results to load." )
			.set_is_required()
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 1 )
			.set_maximum_multiplicity( 100 )
		;
		
		options.declare_group( "SNP exclusion options" ) ;
		options[ "-excl-snpids" ]
			.set_description( "Exclude all SNPs whose SNPID is in the given file(s) from the analysis.")
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-excl-rsids" ]
			.set_description( "Exclude all SNPs whose RSID is in the given file(s) from the analysis.")
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-incl-snpids" ]
			.set_description( "Exclude all SNPs whose SNPID is not in the given file(s) from the analysis.")
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-incl-rsids" ]
			.set_description( "Exclude all SNPs whose RSID is not in the given file(s) from the analysis.")
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-excl-snps-matching" ]
			.set_description( "Filter out snps whose rsid or SNPID matches the given value. "
				"The value should be a string which can contain a % wildcard character (which matches any substring). "
				"Optionally, prefix the argument with snpid~ or rsid~ to only match against the SNPID or rsid fields." )
			.set_takes_single_value() ;
		options[ "-incl-snps-matching" ]
			.set_description( "Filter out snps whose rsid or SNPID does not match the given value. "
				"The value should be a string which can contain a % wildcard character (which matches any substring). "
				"Optionally, prefix the argument with snpid~ or rsid~ to only match against the SNPID or rsid fields." )
			.set_takes_single_value() ;
		options[ "-incl-range" ]
			.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to operate on. "
				"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
				"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
				"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
			.set_takes_single_value() ;
		options[ "-excl-range" ]
			.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to exclude from operation. "
				"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
				"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
				"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
			.set_takes_single_value() ;

		options[ "-min-info" ]
			.set_description( "Treat SNPs with info less than the given threshhold as missing." )
			.set_takes_values( 1 ) ;

		options[ "-min-maf" ]
			.set_description( "Treat SNPs with maf (in controls) less than the given threshhold as missing." )
			.set_takes_values( 1 ) ;
		
		options[ "-o" ]
			.set_description( "Specify the path to the output file." )
			.set_is_required()
			.set_takes_single_value() ;
		
		options[ "-analysis-name" ]
			.set_description( "Specify a name to label results from this analysis with" )
			.set_takes_single_value()
			.set_default_value( "bingwa analysis, started " + appcontext::get_current_time_as_string() ) ;
		
		options[ "-snp-match-fields" ]
			.set_description( "Use this option to specify a comma-separated list of SNP-identifying fields that should be used "
				"to match SNPs between cohorts.  Possible fields "
				"are \"position\", \"alleles\", \"rsid\", or \"snpid\"; you must always specify \"position\" as the first entry." )
			.set_takes_single_value()
			.set_default_value( "position,alleles" ) ;

		options[ "-log" ]
			.set_description( "Specify the path of a log file; all screen output will be copied to the file." )
			.set_takes_single_value() ;

		options.declare_group( "Analysis options" ) ;
		
		options[ "-prior-correlation" ]
			.set_description( "Specify a correlation matrix of the form\n"
				"   [ 1 r r  .. ]\n"
				"   [ r 1 r  .. ]\n"
				"   [ r r 1  .. ]\n"
				"   [       .   ]\n"
				"   [         . ]\n"
				" where r denotes the between-cohort correlation."
			)
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 0 )
			.set_maximum_multiplicity( 1 ) ; // run up to 100 models

		options[ "-prior-correlation-matrix" ]
			.set_description( "Specify the upper triangle of the prior correlation matrix for bayesian analysis.  For example, "
				"the value \"1,0.5,1\" specifies the matrix\n"
				"   [ 1    0.5 ]\n"
				"   [ 0.5  1   ]."
			)
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 0 )
			.set_maximum_multiplicity( 1 ) ; // run up to 100 models

		options[ "-prior-sd" ]
			.set_description( "Specify the prior standard deviation for bayesian analysis." )
			.set_takes_single_value() ;
		
		options.option_implies_option( "-prior-correlation", "-prior-sd" ) ;
		options.option_implies_option( "-prior-correlation-matrix", "-prior-sd" ) ;
		options.option_excludes_option( "-prior-correlation", "-prior-correlation-matrix" ) ;
	}
} ;


namespace impl {
	void assign_counts( Eigen::MatrixXd* result, Eigen::MatrixXd::Index const col, genfile::VariantEntry const& AA, genfile::VariantEntry const& AB, genfile::VariantEntry const& BB, genfile::VariantEntry const& missing ) {
		assert( result ) ;
		assert( result->rows() >= 4 ) ;
		assert( result->cols() > col ) ;
		double const NA =  std::numeric_limits< double >::quiet_NaN() ;
		(*result)( 0, col ) = AA.is_missing() ? NA : AA.as< double >() ;
		(*result)( 1, col ) = AB.is_missing() ? NA : AB.as< double >() ;
		(*result)( 2, col ) = BB.is_missing() ? NA : BB.as< double >() ;
		(*result)( 3, col ) = missing.is_missing() ? NA : missing.as< double >() ;
	}
	
	void assign_vector_elt( Eigen::VectorXd* result, Eigen::VectorXd* non_missingness, Eigen::VectorXd::Index const elt, genfile::VariantEntry const& value ) {
		assert( result ) ;
		assert( non_missingness ) ;
		assert( result->size() == non_missingness->size() ) ;
		assert( result->size() > elt ) ;
		double const NA =  std::numeric_limits< double >::quiet_NaN() ;
		if( value.is_missing() ) {
			(*non_missingness)( elt ) = 0.0 ;
			(*result)( elt ) = NA ;
		}
		else {
			(*result)( elt ) = value.as< double >() ;
		}
	}
}

struct AmetProcessor: public boost::noncopyable
{
	typedef std::auto_ptr< AmetProcessor > UniquePtr ;
	static UniquePtr create( genfile::SNPIdentifyingData2::CompareFields const& compare_fields ) {
		return UniquePtr( new AmetProcessor( compare_fields ) ) ;
	}
	
	AmetProcessor( genfile::SNPIdentifyingData2::CompareFields const& compare_fields ):
		m_snps( compare_fields )
	{
	}
	
	void add_cohort( std::string const& name, FrequentistGenomeWideAssociationResults::UniquePtr results ) {
		assert( results.get() ) ;
		m_cohort_names.push_back( name ) ;
		m_cohorts.push_back( results ) ;
	}
	
	void summarise( appcontext::UIContext& ui_context ) {
		using genfile::string_utils::to_string ;
		
		ui_context.logger() << "================================================\n" ;
		ui_context.logger() << "Cohort summary:\n" ;
		for( std::size_t i = 0; i < m_cohorts.size(); ++i ) {
			ui_context.logger() << " - Cohort " << (i+1) << ":\n" ;
			ui_context.logger() << m_cohorts[i].get_summary( "   - " ) ;
			ui_context.logger() << "\n   - First few SNPs are:\n" ;
			for( std::size_t snp_i = 0; snp_i < std::min( std::size_t( 5 ), m_cohorts[i].get_number_of_SNPs() ); ++snp_i ) {
				double info ;
				double maf ;
				m_cohorts[i].get_info( snp_i, &info ) ;
				m_cohorts[i].get_maf_in_controls( snp_i, &maf ) ;
				ui_context.logger() << "     " << m_cohorts[i].get_SNP( snp_i ) << " (maf = " << maf << ", info = " << info << ")\n";
			}
		}
		
		ui_context.logger() << "\n================================================\n" ;
		ui_context.logger() << "SNP Categories:\n" ;
		ui_context.logger() << "  " ;
		for( std::size_t i = 0; i < m_cohorts.size(); ++i ) {
			ui_context.logger() << std::setw( 12 ) << ( "in scan " + to_string( i+1 )) ;
		}
		ui_context.logger() << "\n" ;
		for( CategoryCounts::const_iterator i = m_category_counts.begin(); i != m_category_counts.end(); ++i ) {
			ui_context.logger() << "  " ;
			for( std::size_t j = 0; j < m_cohorts.size(); ++j ) {	
				ui_context.logger() << std::setw(12) << ( i->first[j] ? '*' : ' ' ) ;
			}
			ui_context.logger() << ": " << i->second << "\n" ;
		}
		ui_context.logger() << " TOTAL: " << m_snps.size() << "\n" ;
		ui_context.logger() << "================================================\n" ;
	}
	
	void add_computation( std::string const& name, AmetComputation::UniquePtr computation ) {
		m_computations.push_back( computation ) ;
	}

	void setup( appcontext::UIContext& ui_context ) {
		unsafe_setup( ui_context ) ;
	}

	void process( appcontext::UIContext& ui_context ) {
		unsafe_process( ui_context ) ;
	}

	typedef boost::signals2::signal<
		void (
			std::size_t index,
			genfile::SNPIdentifyingData2 const& snp,
			std::string const& computation_name,
			std::string const& value_name,
			genfile::VariantEntry const& value
		)
	> ResultSignal ;

	void send_results_to( ResultSignal::slot_type callback ) {
		m_result_signal.connect( callback ) ;
	}
	
private:
	std::vector< std::string > m_cohort_names ;
	boost::ptr_vector< FrequentistGenomeWideAssociationResults > m_cohorts ;
	boost::ptr_vector< AmetComputation > m_computations ;
	typedef std::map< genfile::SNPIdentifyingData2, std::vector< boost::optional< std::size_t > >, genfile::SNPIdentifyingData2::CompareFields > SnpMap ;
	SnpMap m_snps ;
	typedef std::map< std::vector< bool >, std::size_t > CategoryCounts ;
	CategoryCounts m_category_counts ;
	
	ResultSignal m_result_signal ;

	struct DataGetter: public AmetComputation::DataGetter {
		DataGetter(
			boost::ptr_vector< FrequentistGenomeWideAssociationResults > const& cohorts,
			std::vector< boost::optional< std::size_t > >const& indices
		):
			m_cohorts( cohorts ),
			m_indices( indices )
		{}
		
		std::size_t get_number_of_cohorts() const { return m_cohorts.size() ; }
		
		void get_counts( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				return m_cohorts[i].get_counts( *m_indices[i], result ) ;
			}
		}
		void get_betas( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_betas( *m_indices[i], result ) ;
			}
		}
		void get_ses( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_ses( *m_indices[i], result ) ;
			}
		}
		void get_pvalue( std::size_t i, double* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_pvalue( *m_indices[i], result ) ;
			}
		}
		void get_info( std::size_t i, double* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_info( *m_indices[i], result ) ;
			}
		}
		bool is_non_missing( std::size_t i ) const {
			return( m_indices[i] ) ;
		}

		private:
			boost::ptr_vector< FrequentistGenomeWideAssociationResults > const& m_cohorts ;
			std::vector< boost::optional< std::size_t > > const& m_indices ;
		
	} ;
private:
	void unsafe_setup( appcontext::UIContext& ui_context ) {
		link_data( ui_context ) ;
		categorise_snps() ;
	}
	
	void add_SNP_callback( std::size_t cohort_i, std::size_t snp_i, genfile::SNPIdentifyingData2 const& snp ) {
		// Find the SNP that matches the given one (if it exists)
		std::pair< SnpMap::iterator, SnpMap::iterator > range = m_snps.equal_range( snp ) ;
		if( range.second == range.first ) {
			// no match, so add this SNP.
			std::vector< boost::optional< std::size_t > >& snp_indices = m_snps[ snp ] ;
			snp_indices.resize( m_cohorts.size() ) ;
			snp_indices[ cohort_i ] = snp_i ;
		}
		else {
			// There is a match.  Combine the rsid, snpid, and alleles.
			genfile::SNPIdentifyingData2 stored_snp = range.first->first ;
			std::vector< boost::optional< std::size_t > > snp_indices = range.first->second ;
			m_snps.erase( range.first ) ;
			if( stored_snp.get_rsid() != snp.get_rsid() ) {
				stored_snp.set_rsid( std::string( stored_snp.get_rsid() ) + "," + std::string( snp.get_rsid() ) ) ;
			}
			if( stored_snp.get_SNPID() != snp.get_SNPID() ) {
				stored_snp.set_SNPID( std::string( stored_snp.get_SNPID() ) + "," + std::string( snp.get_SNPID() ) ) ;
			}
			if( stored_snp.get_first_allele() != snp.get_first_allele() ) {
				stored_snp.set_first_allele( std::string( stored_snp.get_first_allele() ) + "," + std::string( snp.get_first_allele() ) ) ;
			}
			if( stored_snp.get_second_allele() != snp.get_second_allele() ) {
				stored_snp.set_second_allele( std::string( stored_snp.get_second_allele() ) + "," + std::string( snp.get_second_allele() ) ) ;
			}
			snp_indices[ cohort_i ] = snp_i ;
			m_snps[ stored_snp ] = snp_indices ;
		}
	}

	void link_data( appcontext::UIContext& ui_context ) {
		appcontext::UIContext::ProgressContext progress_context = ui_context.get_progress_context( "Linking SNPs" ) ;
		progress_context( 0, m_cohorts.size() ) ;
		for( std::size_t cohort_i = 0; cohort_i < m_cohorts.size(); ++cohort_i ) {
			m_cohorts[ cohort_i].get_SNPs(
				boost::bind(
					&AmetProcessor::add_SNP_callback,
					this,
					cohort_i,
					_1,
					_2
				)
			) ;
			progress_context( cohort_i+1, m_cohorts.size() ) ;
		}
	}
	
	void categorise_snps() {
		SnpMap::const_iterator snp_i = m_snps.begin(), end_i = m_snps.end() ;
		for( ; snp_i != end_i; ++snp_i ) {
			std::vector< boost::optional< std::size_t > > const& indices = snp_i->second ;
			std::vector< bool > indicator( indices.size(), false ) ;
			for( std::size_t i = 0; i < indices.size(); ++i ) {
				indicator[i] = bool( indices[i] ) ;
			}
			++m_category_counts[ indicator ] ;
		}
	}

	void unsafe_process( appcontext::UIContext& ui_context ) {
		Eigen::MatrixXd cohort_counts( 4, m_cohorts.size() ) ;
		Eigen::MatrixXd cohort_betas( 2, m_cohorts.size() ) ;
		Eigen::MatrixXd cohort_ses( 2, m_cohorts.size() ) ;
		Eigen::VectorXd non_missingness( m_cohorts.size() ) ;

		{
			appcontext::UIContext::ProgressContext progress_context = ui_context.get_progress_context( "Storing variants" ) ;
		
			SnpMap::const_iterator snp_i = m_snps.begin() ;
			SnpMap::const_iterator const end_i = m_snps.end() ;
			for( std::size_t snp_index = 0; snp_i != end_i; ++snp_i, ++snp_index ) {
				m_result_signal( snp_index, snp_i->first, "AmetProcessor", "index", genfile::VariantEntry::Integer( snp_index ) ) ;
				progress_context( snp_index + 1, m_snps.size() ) ;
			}
		}
		{
			appcontext::UIContext::ProgressContext progress_context = ui_context.get_progress_context( "Storing meta-analysis results" ) ;
			SnpMap::const_iterator snp_i = m_snps.begin() ;
			SnpMap::const_iterator const end_i = m_snps.end() ;
			for( std::size_t snp_index = 0; snp_i != end_i; ++snp_i, ++snp_index ) {
				std::vector< boost::optional< std::size_t > > const& indices = snp_i->second ;
				DataGetter data_getter( m_cohorts, indices ) ;
				for( std::size_t i = 0; i < m_computations.size(); ++i ) {
					m_computations[i](
						snp_i->first,
						data_getter,
						boost::bind(
							boost::ref( m_result_signal ),
							snp_index,
							snp_i->first,
							m_computations[i].get_spec(),
							_1, _2
						)
					) ;
				}

				progress_context( snp_index + 1, m_snps.size() ) ;
			}
		}
	}
} ;


struct AmetApplication: public appcontext::ApplicationContext {
public:
	AmetApplication( int argc, char **argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version,
			std::auto_ptr< appcontext::OptionProcessor >( new AmetOptions ),
			argc,
			argv,
			"-log"
		),
		m_processor( AmetProcessor::create( genfile::SNPIdentifyingData2::CompareFields( options().get< std::string > ( "-snp-match-fields" )) ) )
	{}
	
	void run() {
		load_data() ;

		m_processor->setup( get_ui_context() ) ;
		m_processor->summarise( get_ui_context() ) ;
		
		snp_summary_component::DBOutputter::SharedPtr outputter = snp_summary_component::DBOutputter::create_shared(
			options().get< std::string >( "-o" ),
			options().get< std::string >( "-analysis-name" ),
			options().get_values_as_map()
		) ;

		m_processor->send_results_to(
			boost::bind(
				&snp_summary_component::DBOutputter::operator(),
				outputter,
				_1, _2, _3, _4, _5
			)
		) ;

		if( options().check( "-snptest" ) && options().get_values< std::string >( "-snptest" ).size() > 1 ) {
			m_processor->add_computation(
				"PerCohortValueReporter",
				AmetComputation::create( "PerCohortValueReporter", options() )
			) ;
			m_processor->add_computation(
				"FixedEffectFrequentistMetaAnalysis",
				AmetComputation::create( "FixedEffectFrequentistMetaAnalysis", options() )
			) ;
			if( options().check( "-prior-correlation" ) || options().check( "-prior-correlation-matrix" )) {
				std::map< std::string, Eigen::MatrixXd > const priors = get_priors( options() ) ;
				assert( priors.size() > 0 ) ;
				std::map< std::string, Eigen::MatrixXd >::const_iterator i = priors.begin() ;
				std::map< std::string, Eigen::MatrixXd >::const_iterator const end_i = priors.end() ;
				for( ; i != end_i; ++i ) {
					AmetComputation::UniquePtr computation(
						new ApproximateBayesianMetaAnalysis(
							i->first,
							i->second
						)
					) ;
					m_processor->add_computation(
						"ApproximateBayesianMetaAnalysis",
						computation
					) ;
				}
			}
		}
		m_processor->process( get_ui_context() ) ;
	}
	
	std::map< std::string, Eigen::MatrixXd > get_priors( appcontext::OptionProcessor const& options ) {
		std::map< std::string, Eigen::MatrixXd > result ;
		using genfile::string_utils::to_string ;
		using genfile::string_utils::to_repr ;
		using genfile::string_utils::split_and_strip_discarding_empty_entries ;
		int const N = options.get_values< std::string >( "-snptest" ).size() ;
		std::vector< std::string > sds = split_and_strip_discarding_empty_entries( options.get< std::string >( "-prior-sd" ), ",", " \t" ) ;
		
		if( options.check( "-prior-correlation" ) ) {
			std::vector< std::string > const rhos = split_and_strip_discarding_empty_entries( options.get< std::string >( "-prior-correlation" ), ",", " \t" ) ;
			if( sds.size() != rhos.size() ) {
				throw genfile::BadArgumentError(
					"AmetProcessor::get_priors()",
					"-prior-correlation \"" + options.get< std::string >( "-prior-correlation" ) + "\""
				) ;
			}
			for( std::size_t i = 0; i < sds.size(); ++i ) {
				double const rho = to_repr< double >( rhos[i] ) ;
				double const sd = to_repr< double >( sds[i] ) ;
				result[ "rho=" + rhos[i] + "/sd=" + sds[i] ] = get_prior_matrix( N, rho, sd ) ;
			}
		} else if( options.check( "-prior-correlation-matrix" )) {
			std::vector< std::string > const matrix_specs = options.get_values< std::string >( "-prior-correlation-matrix" ) ;
			if( sds.size() != matrix_specs.size() ) {
				throw genfile::BadArgumentError(
					"AmetProcessor::get_priors()",
					"-prior-correlation-matrix \"" + options.get< std::string >( "-matrix-specs" ) + "\""
				) ;
			}
			for( std::size_t i = 0; i < sds.size(); ++i ) {
				double const sd = to_repr< double >( sds[i] ) ;
				result[ "model_" + to_string( i+1 ) + "/sd=" + sds[i] ] = get_prior_matrix( 
					N,
					matrix_specs[i],
					sd
				) ;
			}
		}

		return result ;
	}

	Eigen::MatrixXd get_prior_matrix( int const n, double const rho, double const sd ) const {
		return get_correlation_matrix( n, rho ) * sd * sd ;
	}

	Eigen::MatrixXd get_prior_matrix( int const n, std::string const& matrix_spec, double const sd ) const {
		return parse_correlation_matrix( n, matrix_spec ) * sd * sd ;
	}

	Eigen::MatrixXd get_correlation_matrix( int const n, double rho ) const {
		Eigen::MatrixXd result = Eigen::MatrixXd::Identity( n, n ) ;
		for( int i = 0; i < (n-1); ++i ) {
			for( int j = (i+1); j < n; ++j ) {
				result( i, j ) = result( j, i ) = rho ;
			}
		}
		return result ;
	}

	Eigen::MatrixXd parse_correlation_matrix( int const n, std::string const& matrix_spec ) const {
		Eigen::MatrixXd result = Eigen::MatrixXd::Zero( n, n ) ;
		std::vector< std::string > values = genfile::string_utils::split_and_strip( matrix_spec, " " ) ;
		if( values.size() != (( n * (n+1) ) / 2 ) ) {
			throw genfile::BadArgumentError( "AmetProcessor::parse_correlation_matrix()", "matrix_spec=\"" + matrix_spec + "\"" ) ;
		}

		{
			int index = 0 ;
			for( int i = 0; i < n; ++i ) {
				for( int j = 0; j < n; ++j ) {
					result( i, j )
						= ( j >= i ) ?
							genfile::string_utils::to_repr< double >( values[ index++ ] )
							:
							result( j, i ) ;
				}
			}
		}
		return result ;
	}

	void post_summarise() const {
		std::string const output_file = options().get< std::string >( "-o" ) ;
		std::string const analysis = options().get< std::string >( "-analysis-name" ) ;
		std::string const SQL = "SELECT * FROM SummaryDataView WHERE analysis == '" + analysis + "'" ;
		get_ui_context().logger() << "\n" << globals::program_name << ": Analysis complete.\n" ;
		get_ui_context().logger() << "\n---- VIEWING RESULTS ----\n" ;
		get_ui_context().logger() << "View your results from the command-line using this command\n\n"
			<< "$ sqlite3 -separator $'\\t' -header \"" << output_file << "\" \"" << SQL << "\"\n\n"
			<< "or from R using this snippet:\n\n"
			<< "> library( RSQLite ) ;\n"
			<< "> db = dbConnect( dbDriver( \"SQLite\" ), \"" << output_file << "\" ) ;\n"
			<< "> data = dbGetQuery( db, \"" << SQL << "\")  ;\n\n"
			<< "...enjoy!\n" ;
		get_ui_context().logger() << "\n-------------------------\n" ;
		
	}
	
	typedef
		boost::function< void ( std::size_t, genfile::SNPIdentifyingData2 const&, std::string const&, std::string const&, genfile::VariantEntry const& ) > 
		ResultCallback ;

	void load_data() {
		using genfile::string_utils::to_string ;

		std::vector< std::string > cohort_files = options().get_values< std::string >( "-snptest" ) ;
		for( std::size_t i = 0; i < cohort_files.size(); ++i ) {
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading SNPTEST results \"" + cohort_files[i] + "\"" ) ;
			FrequentistGenomeWideAssociationResults::UniquePtr results ;

			SNPTESTResults::SNPResultCallback snp_callback ;
			if( cohort_files.size() == 1 ) {
				snp_summary_component::DBOutputter::SharedPtr raw_outputter = snp_summary_component::DBOutputter::create_shared(
					options().get< std::string >( "-o" ),
					globals::program_name + " cohort " + to_string( i+1 ) + ": SNPTEST analysis: \"" + cohort_files[i] + "\"",
					snp_summary_component::DBOutputter::Metadata()
				) ;

				snp_callback = boost::bind(
					&snp_summary_component::DBOutputter::operator(),
					raw_outputter,
					_1, _2, _3, _4, _5
				) ;
			}

			SNPExclusionTestConjunction::UniquePtr test( new SNPExclusionTestConjunction() ) ;

			if( options().check_if_option_was_supplied_in_group( "SNP exclusion options" ) ) {
				genfile::CommonSNPFilter::UniquePtr snp_filter = get_snp_exclusion_filter() ;
				if( snp_filter.get() ) {
					test->add_subtest( genfile::SNPIdentifyingDataTest::UniquePtr( snp_filter.release() ) ) ;
				}
				if( options().check( "-min-info" ) || options().check( "-min-maf" ) ) {
					ThreshholdingSummaryStatisticTest::UniquePtr subtest( new ThreshholdingSummaryStatisticTest() ) ;
					if( options().check( "-min-info" ) ) {
						subtest->set_inclusion_bounds( "info", options().get< double >( "-min-info" ), 2 ) ; // info should be no larger than 1, but I guess there may be numerical error.
					}
					if( options().check( "-min-maf" ) ) {
						subtest->set_inclusion_bounds( "controls_maf", options().get< double >( "-min-maf" ), 1 ) ;
					}
					test->add_subtest( SNPExclusionTest::UniquePtr( subtest.release() ) ) ;
				}
			}

			results.reset(
				new SNPTESTResults(
					genfile::wildcard::find_files_by_chromosome( cohort_files[i] ),
					SNPExclusionTest::UniquePtr( test.release() ),
					snp_callback,
					progress_context
				)
			) ;
			
			// summarise right now so as to see memory used.
			get_ui_context().logger() << "Cohort " << (i+1) << " summary: " << results->get_summary() << ".\n" ;
			
			m_processor->add_cohort( "cohort_" + to_string( i+1 ), results ) ;
		}
	}
	
	genfile::CommonSNPFilter::UniquePtr get_snp_exclusion_filter() const {
		genfile::CommonSNPFilter::UniquePtr snp_filter ;

		if( options().check_if_option_was_supplied_in_group( "SNP exclusion options" )) {
			snp_filter.reset( new genfile::CommonSNPFilter ) ;

			if( options().check_if_option_was_supplied( "-excl-snpids" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-snpids" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::SNPIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-snpids" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-snpids" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::SNPIDs
					) ;
				}
			}


			if( options().check_if_option_was_supplied( "-excl-rsids" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-rsids" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-rsids" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-rsids" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-snps-matching" )) {
				std::string it = options().get< std::string > ( "-excl-snps-matching" ) ;
				std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( it, ",", " " ) ;
				BOOST_FOREACH( std::string const& spec, specs ) {
					snp_filter->exclude_snps_matching(
						spec
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-snps-matching" )) {
				std::string it = options().get< std::string > ( "-incl-snps-matching" ) ;
				std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( it, ",", " " ) ;
				BOOST_FOREACH( std::string const& spec, specs ) {
					snp_filter->include_snps_matching(
						spec
					) ;
				}
			}
			
			if( options().check_if_option_was_supplied( "-incl-range" )) {
				std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( options().get_value< std::string >( "-incl-range" ), ",", " \t" ) ;
				for ( std::size_t i = 0; i < specs.size(); ++i ) {
					snp_filter->include_snps_in_range(
						genfile::GenomePositionRange::parse( specs[i] )
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-range" )) {
				std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( options().get_value< std::string >( "-excl-range" ), ",", " \t" ) ;
				for ( std::size_t i = 0; i < specs.size(); ++i ) {
					snp_filter->exclude_snps_in_range(
						genfile::GenomePositionRange::parse( specs[i] )
					) ;
				}
			}
		}
		return snp_filter ;
	}
	
private:
	AmetProcessor::UniquePtr m_processor ;
} ;

int main( int argc, char **argv ) {
	try {
		std::cerr << "sizeof( std::size_t ) = " << sizeof( std::size_t ) << ".\n" ;
		std::cerr << "sizeof( std::string ) = " << sizeof( std::string ) << ".\n" ;
		std::cerr << "sizeof( genfile:string_utils::slice ) = " << sizeof( genfile::string_utils::slice ) << ".\n" ;
		std::cerr << "sizeof( genfile::SNPIdentifyingData ) = " << sizeof( genfile::SNPIdentifyingData ) << ".\n" ;
		std::cerr << "sizeof( genfile::SNPIdentifyingData2 ) = " << sizeof( genfile::SNPIdentifyingData2 ) << ".\n" ;
		AmetApplication app( argc, argv ) ;	
		app.run() ;
		app.post_summarise() ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
