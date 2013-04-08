
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
#include "qcdb/Storage.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"
#include "qcdb/FlatFileOutputter.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"

namespace globals {
	std::string const program_name = "bingwa" ;
	std::string const program_version = "0.2" ;
}

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
			if( variable != "maf" && variable != "info" ) {
				throw genfile::BadArgumentError( "ThreshholdingSummaryStatisticTest::set_inclusion_bounds()", "variable=\"" + variable + "\"" ) ;
			}
			m_bounds[ variable ] = std::make_pair( lower_bound, upper_bound ) ;
		}
		
		bool operator()( genfile::SNPIdentifyingData const& snp, double const maf, double const info ) const {
			VariableBounds::const_iterator i = m_bounds.begin(), end_i = m_bounds.end() ;
			double value ;
			for( ; i != end_i; ++i ) {
				if( i->first == "maf" ) {
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

struct FrequentistGenomeWideAssociationResults: public boost::noncopyable {
public:
	typedef std::auto_ptr< FrequentistGenomeWideAssociationResults > UniquePtr ;
	typedef boost::function< void ( std::size_t i, genfile::SNPIdentifyingData2 const& snp ) > GetSNPCallback ;
	typedef
		boost::function< void ( genfile::SNPIdentifyingData2 const&, std::string const&, genfile::VariantEntry const& ) > 
		SNPResultCallback ;
	typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;

	static UniquePtr create(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		std::vector< std::string > const& columns,
		SNPExclusionTest::UniquePtr test,
		SNPResultCallback callback = SNPResultCallback(),
		ProgressCallback progress_callback = ProgressCallback()
	) ;
public:
	virtual ~FrequentistGenomeWideAssociationResults() {}
	
	virtual void add_variable( std::string const& variable ) = 0 ;
	
	virtual std::size_t get_number_of_SNPs() const = 0 ;
	virtual genfile::SNPIdentifyingData2 const& get_SNP( std::size_t snp_i ) const = 0 ;
	virtual void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ;
	virtual void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ; 
	virtual void get_pvalue( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ;
	virtual void get_info( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_frequency( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_variable( std::size_t snp_i, std::string const& variable, double* result ) const = 0 ;
	
	virtual std::string get_summary( std::string const& prefix = "", std::size_t target_column = 80 ) const = 0 ;

	// Legacy function: this should probably be removed.
	void get_SNPs( GetSNPCallback callback ) const {
		for( std::size_t i = 0; i < get_number_of_SNPs(); ++i ) {
			callback( i, get_SNP( i ) ) ;
		}
	}
} ;

/*
* Class FlatFileScanResults.
* This is a utility base class for classes reading scan results flat files.
* It provides:
* - storage for values read.
* - 
*
* - and implements the FrequentistGenomeWideAssociationResults interface for retrieving tresults..
*/
struct FlatFileScanResults:  public FrequentistGenomeWideAssociationResults {
public:
	typedef std::auto_ptr< FlatFileScanResults > UniquePtr ;
public:
	void add_data(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) {
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
	void get_frequency( std::size_t snp_i, double* result ) const {
		*result = ( 2.0 * m_sample_counts( snp_i, 2 ) + m_sample_counts( snp_i, 1 ) ) / ( 2.0 * m_sample_counts.row( snp_i ).sum() ) ;
	}
	void get_variable( std::size_t snp_i, std::string const& variable, double* value ) const {
		std::map< std::string, std::vector< double > >::const_iterator where = m_extra_variables.find( variable ) ;
		assert( where != m_extra_variables.end() ) ;
		*value = where->second[ snp_i ] ;
	}

	std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
		using genfile::string_utils::to_string ;
		
		// estimate memory used in SNPs.
		unsigned long mem_used = 0 ;
		for( std::size_t i = 0; i < m_snps.size(); ++i ) {
			mem_used += m_snps[i].get_estimated_bytes_used() ;
		}
		
		std::string result = "scan results object ("
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
						+ m_maf.size()
						+ m_sample_counts.size()
					) * sizeof( float )
					+ mem_used
				)
				/ 1000000.0
			) + "Mb in use.)" ;
		return result ;
	}

protected:
	std::string const m_missing_value ;
	SNPExclusionTest::UniquePtr m_exclusion_test ;
	typedef boost::bimap< std::string, std::size_t > ColumnMap ;
	
	std::vector< genfile::SNPIdentifyingData2 > m_snps ;
	Eigen::MatrixXf m_betas ;
	Eigen::MatrixXf m_ses ;
	Eigen::VectorXf m_pvalues ;
	Eigen::VectorXf m_info ;
	Eigen::VectorXf m_maf ;
	Eigen::MatrixXf m_sample_counts ;
	typedef std::map< std::string, std::vector< double > > ExtraVariables ;
	ExtraVariables m_extra_variables ;
	
	FlatFileScanResults():
		m_missing_value( "NA" )
	{}
	
	virtual std::set< std::string > get_desired_columns() const = 0 ;
	virtual bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const = 0 ;
	virtual bool check_if_snp_accepted( std::size_t snp_index ) const = 0 ;
	virtual void store_value( int snp_index, std::string const& variable, double value ) = 0 ;
	
private:

	void setup(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) {
		progress_callback( 0, 100 ) ;
		setup( statfile::BuiltInTypeStatSource::open( filenames ), callback, progress_callback ) ;
	}

	void setup(
		statfile::BuiltInTypeStatSource::UniquePtr source,
		SNPResultCallback callback,
		ProgressCallback progress_callback
	) {
		ColumnMap column_map = get_columns_to_store( *source ) ;
		
		int degrees_of_freedom = ( column_map.left.find( "_beta_2" ) == column_map.left.end() ) ? 1 : 2 ;
		
		double const NA = std::numeric_limits< double >::quiet_NaN() ;
		
		std::size_t snp_index = m_snps.size() ;

		genfile::SNPIdentifyingData snp ;
		Eigen::VectorXd betas ;
		Eigen::VectorXd ses ;
		
		for( ; read_snp( *source, snp ); (*source) >> statfile::ignore_all() ) {
			
			if( snp_index >= m_snps.size() ) {
				resize_storage( snp_index + 100000, degrees_of_freedom ) ;
			}
			
			// Deal with strange non-ids.  This isn't a general solution but 
			if(
				snp.get_SNPID() == "---" // IMPUTE2
				|| snp.get_SNPID() == "?" // QCTOOL under some usages
				|| snp.get_SNPID() == "NA" // dunno.
			) {
				snp.set_SNPID( "" ) ;
			}

			// read data columns

			ColumnMap::right_const_iterator i = column_map.right.begin(), end_i = column_map.right.end() ;
			std::string value ;
			for( ; i != end_i; ++i ) {
				(*source)
					>> statfile::ignore( i->first - source->current_column() )
					>> value ;

				store_value(
					snp_index,
					i->second,
					( value == m_missing_value ) ? NA : genfile::string_utils::to_repr< double >( value )
				) ;
			}

			m_snps[ snp_index ] = snp ;
			if( check_if_snp_accepted( snp_index ) ) {
				++snp_index ;
			}

			if( progress_callback ) {
				progress_callback( double( source->number_of_rows_read() + 1 ), source->number_of_rows() ) ;
			}
		}
		
		// Now deallocate any unused memory.
		m_snps.resize( snp_index ) ;
		resize_storage( snp_index, degrees_of_freedom ) ;
	}
	
	ColumnMap get_columns_to_store(
		statfile::BuiltInTypeStatSource const& source
	) {
		using genfile::string_utils::to_string ;
		std::set< std::string > desired_columns = get_desired_columns() ;

		ColumnMap result ;
		for( std::size_t i = 0; i < source.number_of_columns(); ++i ) {
			std::string name = source.name_of_column( i ) ;
			for( std::set< std::string >::iterator j = desired_columns.begin(); j != desired_columns.end(); ++j ) {
				assert( j->size() > 0 ) ;
				if(
					(
						j->at(0) == '*'
						&& name.size() >= ( j->size() - 1 )
						&& name.compare( name.size() - j->size() + 1, j->size() - 1, j->substr( 1, j->size() ) ) == 0
					) || (
						j->at( j->size() - 1 ) == '*'
						&& name.size() >= ( j->size() - 1 )
						&& name.compare( 0, j->size() - 1, j->substr( 0, j->size() - 1 ) ) == 0
					) || (
						j->at(0) != '*' && j->at( j->size() - 1 ) != '*' && name == *j
					)
				) {
					std::pair< ColumnMap::iterator, bool > insertion = result.insert( ColumnMap::value_type( *j, i )) ;
					if( !insertion.second ) {
						throw genfile::OperationFailedError(
							"SNPTESTResults::get_columns_to_store()",
							"m_column_map",
							"Insertion of value for column \"" + name + "\" (matching \"" + to_string( *j ) + "\")."
						) ;
					}
				}
			}
		}
		
		std::set< std::string > required_columns = desired_columns ;
		required_columns.erase( "*_beta_2" ) ;
		required_columns.erase( "*_se_2" ) ;
		
		for( std::set< std::string >::const_iterator i = required_columns.begin(); i != required_columns.end(); ++i ) {
			if( result.left.find( *i ) == result.left.end() ) {
				throw genfile::BadArgumentError(
					"FlatFileScanResults::get_columns_to_store()",
					"required column=\"" + *i + "\"",
					"Could not find matching column in source \"" + source.get_source_spec() + "\""
				) ;
			}
		}

		return result ;
	}

	void resize_storage( Eigen::MatrixXf::Index const N_snps, Eigen::MatrixXf::Index const degrees_of_freedom ) {
		using std::min ;
		int const current_N = min( N_snps, Eigen::MatrixXf::Index( m_snps.size() ) ) ;
		{
			// free any unused memory for SNPs.
			std::vector< genfile::SNPIdentifyingData2 > snps( N_snps ) ;
			std::copy( m_snps.begin(), m_snps.end(), snps.begin() ) ;
			m_snps.swap( snps ) ;
		}
		{
			Eigen::MatrixXf betas = Eigen::MatrixXf::Zero( N_snps, degrees_of_freedom ) ;
			if( m_betas.rows() > 0 ) {
				betas.block( 0, 0, current_N, degrees_of_freedom ) = m_betas.block( 0, 0, current_N, degrees_of_freedom ) ;
			}
			m_betas.swap( betas ) ;
		}
		{
			Eigen::MatrixXf ses = Eigen::MatrixXf::Zero( N_snps, degrees_of_freedom )  ;
			if( m_ses.rows() > 0 ) {
				ses.block( 0, 0, current_N, degrees_of_freedom ) = m_ses.block( 0, 0, current_N, degrees_of_freedom ) ;
			}
			m_ses.swap( ses ) ;
		}
		{
			Eigen::VectorXf pvalues = Eigen::VectorXf::Zero( N_snps ) ;
			pvalues.head( current_N ) = m_pvalues.head( current_N ) ;
			m_pvalues.swap( pvalues ) ;
		}
		{
			Eigen::VectorXf info = Eigen::VectorXf::Zero( N_snps ) ;
			info.head( current_N ) = m_info.head( current_N ) ;
			m_info.swap( info ) ;
		}
		{
			Eigen::VectorXf maf = Eigen::VectorXf::Zero( N_snps ) ;
			maf.head( current_N ) = m_maf.head( current_N ) ;
			m_maf.swap( maf ) ;
		}
		{
			Eigen::MatrixXf sample_counts = Eigen::MatrixXf::Zero( N_snps, 4 ) ;
			if( m_sample_counts.rows() > 0 ) {
				sample_counts.block( 0, 0, current_N, 4 ) = m_sample_counts.block( 0, 0, current_N, 4 ) ;
			}
			m_sample_counts.swap( sample_counts ) ;
		}
		{
			for( ExtraVariables::iterator i = m_extra_variables.begin(); i != m_extra_variables.end(); ++i ) {
				i->second.resize( N_snps ) ; // std::vector resize does not lose data.
				std::vector< double > v = i->second ;
				i->second.swap( v ) ;
			}
		}
	}
	
	void free_unused_memory() {
		std::size_t const N_snps = m_snps.size() ;
		{
			std::vector< genfile::SNPIdentifyingData2 > snps( m_snps.begin(), m_snps.end() ) ;
			m_snps.swap( snps ) ;
		}
		{
			Eigen::MatrixXf betas = m_betas.block( 0, 0, N_snps, m_betas.cols() ) ;
			m_betas.swap( betas ) ;
		}
		{
			Eigen::MatrixXf ses = m_ses.block( 0, 0, N_snps, m_betas.cols() ) ; ;
			m_ses.swap( ses ) ;
		}
		{
			Eigen::VectorXf pvalues = m_pvalues.head( N_snps ) ;
			m_pvalues.swap( pvalues ) ;
		}
		{
			Eigen::VectorXf info = m_info.head( N_snps ) ;
			m_info.swap( info ) ;
		}
		{
			Eigen::VectorXf maf = m_maf.head( N_snps ) ;
			m_maf.swap( maf ) ;
		}
		{
			Eigen::MatrixXf sample_counts = m_sample_counts.block( 0, 0, N_snps, m_sample_counts.cols() ) ;
			m_sample_counts.swap( sample_counts ) ;
		}
		{
			for( ExtraVariables::iterator i = m_extra_variables.begin(); i != m_extra_variables.end(); ++i ) {
				i->second.resize( N_snps ) ;
				std::vector< double > v = i->second ;
				i->second.swap( v ) ;
			}
		}
	}
} ;

struct SNPTESTResults: public FlatFileScanResults {
	SNPTESTResults(
		SNPExclusionTest::UniquePtr test
	):
	 	m_exclusion_test( test )
	{}

	void add_variable( std::string const& variable ) {
		m_variables.insert( variable ) ;
		/* Make sure we prepare storage. */
		m_extra_variables[ variable ] ;
	}

	
	std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
		return prefix + "SNPTEST " + FlatFileScanResults::get_summary( "", target_column ) ;
	}

private:
	SNPExclusionTest::UniquePtr m_exclusion_test ;
	std::set< std::string > m_variables ;
	
	std::set< std::string > get_desired_columns() const {
		std::set< std::string > desired_columns ;
		desired_columns.insert( "*_beta_1" ) ;
		desired_columns.insert( "*_beta_2" ) ;
		desired_columns.insert( "*_se_1" ) ;
		desired_columns.insert( "*_se_2" ) ;
		desired_columns.insert( "*_pvalue" ) ;
		desired_columns.insert( "info" ) ;
		desired_columns.insert( "all_maf" ) ;
		desired_columns.insert( "all_AA" ) ;
		desired_columns.insert( "all_AB" ) ;
		desired_columns.insert( "all_BB" ) ;
		desired_columns.insert( "all_NULL" ) ;
		desired_columns.insert( m_variables.begin(), m_variables.end() ) ;

		return desired_columns ;
	}

	bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const {
		return( source >> snp.SNPID() >> snp.rsid() >> snp.position().chromosome() >> snp.position().position() >> snp.first_allele() >> snp.second_allele() ) ;
	}
	
	bool check_if_snp_accepted( std::size_t snp_i ) const {
		return
			! m_exclusion_test.get()
			|| m_exclusion_test->operator()( m_snps[ snp_i ], m_maf( snp_i ), m_info( snp_i ) )
		;
	}
	
	void store_value(
		int snp_index,
		std::string const& variable,
		double value
	) {
		using genfile::string_utils::to_repr ;
		
		if( variable == "*_beta_1" ) {
			m_betas( snp_index, 0 ) = value ;
		}
		else if( variable == "*_se_1" ) {
			m_ses( snp_index, 0 ) = value ;
		}
		if( variable == "*_beta_2" ) {
			m_betas( snp_index, 1 ) = value ;
		}
		else if( variable == "*_se_2" ) {
			m_ses( snp_index, 1 ) = value ;
		}
		else if( variable == "*_pvalue" ) {
			m_pvalues( snp_index ) = value ;
		}
		else if( variable == "info" ) {
			m_info( snp_index ) = value ;
		}
		else if( variable == "all_maf" ) {
			m_maf( snp_index ) = value ;
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
		else if( m_variables.find( variable ) != m_variables.end() ) {
			m_extra_variables[ variable ][ snp_index ] = value ;
		}
	}
} ;

struct MMMResults: public FlatFileScanResults {
	MMMResults(
		SNPExclusionTest::UniquePtr test
	):
	 	m_exclusion_test( test )
	{}

	void add_variable( std::string const& variable ) {
		m_variables.insert( variable ) ;
	};

	std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
		return prefix + "mmm " + FlatFileScanResults::get_summary( "", target_column ) ;
	}

private:
	SNPExclusionTest::UniquePtr m_exclusion_test ;
	std::set< std::string > m_variables ;
	
	bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const {
		return( source >> snp.position().chromosome() >> snp.SNPID() >> snp.rsid() >> snp.position().position() >> snp.first_allele() >> snp.second_allele() ) ;
	}

	bool check_if_snp_accepted( std::size_t snp_i ) const {
		return
			!m_exclusion_test.get()
			|| m_exclusion_test->operator()( m_snps[ snp_i ], m_maf( snp_i ), m_info( snp_i ) )
		;
	}
	
	std::set< std::string > get_desired_columns() const {
		std::set< std::string > required_columns ;
		required_columns.insert( "z_est" ) ;
		required_columns.insert( "z_se" ) ;
		required_columns.insert( "pval" ) ;
		required_columns.insert( "var_info_all" ) ;
		required_columns.insert( "freq_1" ) ;
		required_columns.insert( "gen_00" ) ;
		required_columns.insert( "gen_01" ) ;
		required_columns.insert( "gen_11" ) ;
		required_columns.insert( "gen_NULL" ) ;
		required_columns.insert( m_variables.begin(), m_variables.end() ) ;
		return required_columns ;
	}
	
	void store_value(
		int snp_index,
		std::string const& variable,
		double const value
	) {
		using genfile::string_utils::to_repr ;
		
		if( variable == "z_est" ) {
			m_betas( snp_index, 0 ) = value ;
		}
		else if( variable == "z_se" ) {
			m_ses( snp_index, 0 ) = value ;
		}
		else if( variable == "pval" ) {
			m_pvalues( snp_index ) = value ;
		}
		else if( variable == "var_info_all" ) {
			m_info( snp_index ) = value ;
		}
		else if( variable == "freq_1" ) {
			m_maf( snp_index ) = value ;
		}
		else if( variable == "gen_00" ) {
			m_sample_counts( snp_index, 0 ) = value ;
		}
		else if( variable == "gen_01" ) {
			m_sample_counts( snp_index, 1 ) = value ;
		}
		else if( variable == "gen_11" ) {
			m_sample_counts( snp_index, 2 ) = value ;
		}
		else if( variable == "gen_NULL" ) {
			m_sample_counts( snp_index, 3 ) = value ;
		}
	}
} ;

FrequentistGenomeWideAssociationResults::UniquePtr FrequentistGenomeWideAssociationResults::create(
	std::vector< genfile::wildcard::FilenameMatch > const& filenames,
	std::vector< std::string > const& columns,
	SNPExclusionTest::UniquePtr test,
	SNPResultCallback result_callback,
	ProgressCallback progress_callback
) {
	// peek at the file to determine file type
	
	std::string type = "unknown" ;

	if( filenames.size() > 0 ) {
		std::auto_ptr< std::istream > file = genfile::open_text_file_for_input( filenames[0].filename() ) ;
		std::string line ;
		std::getline( *file, line ) ;
		if( line.substr( 0, 32 ) == "chr snp_id1 snp_id2 pos allele_0" ) {
			type = "mmm" ; // Matti's mixed model, http://www.well.ox.ac.uk/~mpirinen/
		}
		else if( line.substr( 0, 40 ) == "id rsid chromosome pos allele_A allele_B" ) {
				type = "snptest" ;
		}
	}

	FlatFileScanResults::UniquePtr result ;
	if( type == "snptest" || type == "unknown" ) {
		result.reset( new SNPTESTResults( test ) ) ;
	}
	else if( type == "mmm" ) {
		result.reset( new MMMResults( test ) ) ;
	}

	BOOST_FOREACH( std::string const& column, columns ) {
		result->add_variable( column ) ;
	}

	result->add_data(
		filenames,
		result_callback,
		progress_callback
	) ;

	return FrequentistGenomeWideAssociationResults::UniquePtr( result.release() ) ;
}

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
		virtual void get_variable( std::string const& variable, std::size_t i, double* result ) const = 0 ;
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
public:
	typedef std::auto_ptr< PerCohortValueReporter > UniquePtr ;
	
public:
	PerCohortValueReporter( std::vector< std::string > const& cohort_names ):
		m_cohort_names( cohort_names )
	{}
	
	void add_variable( std::string const& variable ) {
		m_extra_variables.insert( variable ) ;
	}
	
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
				std::string prefix = m_cohort_names[ i ] + ":" ;
				callback( prefix + "AA", counts(0) ) ;
				callback( prefix + "AB", counts(1) ) ;
				callback( prefix + "BB", counts(2) ) ;
				callback( prefix + "NULL", counts(3) ) ;
				callback( prefix + "B_allele_frequency", ( 2.0 * counts(2) + counts(1) ) / ( 2.0 * counts.head(3).sum() ) ) ;
				
				assert( betas.size() == ses.size() ) ;
				for( int j = 0; j < betas.size(); ++j ) {
					callback( prefix + "beta_" + to_string( j+1 ), betas(j) ) ;
					callback( prefix + "se_" + to_string( j+1 ), ses(j) ) ;
				}
				callback( prefix + "pvalue", pvalue ) ;
				callback( prefix + "info", info ) ;
				{
					double value ;
					BOOST_FOREACH( std::string const& variable, m_extra_variables ) {
						data_getter.get_variable( variable, i, &value ) ;
						callback( prefix + variable, value ) ;
					}
				}
			}
		}
	}
	
	std::string get_spec() const {
		return "PerCohortValueReporter" ;
	}
	
	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}

	private:
		std::vector< std::string > const m_cohort_names ;
		std::set< std::string > m_extra_variables ;
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
		m_sigma( sigma ),
		m_compute_posterior_mean_and_variance( false )
	{
		assert( m_sigma.rows() == m_sigma.cols() ) ;
	}

	void operator()(
		SNPIdentifyingData2 const& snp,
		DataGetter const& data_getter,
		ResultCallback callback
	) ;

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
	bool const m_compute_posterior_mean_and_variance ;

	void compute_bayes_factor( Eigen::MatrixXd const& prior, Eigen::MatrixXd const& V, Eigen::VectorXd const& betas, ResultCallback callback ) const ;

/*
	bayesian_meta_analysis <- function( betas, ses, prior ) {
		V = diag( ses^2 ) ;
		print(V)
		betas = matrix( betas, ncol = 1, nrow = length( betas )) ;
		constant = sqrt( det( V ) / det( V + prior ) ) ;
		exponent = 0.5 * t( betas ) %*% ( solve( V ) - solve( V + prior ) ) %*% betas
		return( constant * exp( exponent ))
	}
*/
} ;

void ApproximateBayesianMetaAnalysis::operator()(
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
		callback( m_prefix + "/bf", genfile::MissingValue() ) ;
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

		compute_bayes_factor( prior, V, betas, callback ) ;
	}
}

namespace impl {
	template< typename M >
	std::string format_matrix( M const& m ) {
		std::ostringstream s ;
		s << "matrix( nrow=" << m.rows() << ", ncol=" << m.cols() << ", data = c(" ;
		for( int i = 0; i < m.rows(); ++i ) {
			if( i > 0 ) {
				s << "," ;
			}
			for( int j = 0; j < m.cols(); ++j ) {
				if( j > 0 ) {
					s << "," ;
				}
				s << m(i,j) ;
			}
		}
		s << "))" ;
		return s.str() ;
	}
}

void ApproximateBayesianMetaAnalysis::compute_bayes_factor( Eigen::MatrixXd const& prior, Eigen::MatrixXd const& V, Eigen::VectorXd const& betas, ResultCallback callback ) const {
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

	if( m_compute_posterior_mean_and_variance ) {
		callback( m_prefix + "/posterior_mean", impl::format_matrix( betas - ( V * V_plus_prior_solver.solve( betas ) ) ) ) ;
		callback( m_prefix + "/posterior_variance", impl::format_matrix( V - V * V_plus_prior_solver.solve( Eigen::MatrixXd::Identity( betas.size(), betas.size() ) ) * V ) ) ;
	}
}

AmetComputation::UniquePtr AmetComputation::create( std::string const& name, appcontext::OptionProcessor const& options ) {
	AmetComputation::UniquePtr result ;
	if( name == "FixedEffectFrequentistMetaAnalysis" ) {
		result.reset( new FixedEffectFrequentistMetaAnalysis() ) ;
	}
	else if( name == "PerCohortValueReporter" ) {
		std::vector< std::string > names ;
		if( options.check( "-cohort-names" )) {
			names = options.get_values< std::string >( "-cohort-names" ) ;
			if( names.size() != options.get_values< std::string >( "-data" ).size() ) {
				throw genfile::BadArgumentError( "AmetComputation::create()", "-cohort-names=\"" + genfile::string_utils::join( names, " " ) + "\"" ) ;
			}
		}
		else {
			names = options.get_values< std::string >( "-data" ) ;
			for( std::size_t i = 0; i < names.size(); ++i ) {
				names[i] = "cohort " + genfile::string_utils::to_string( i + 1 ) + " (file://" + names[i] + ")" ;
			}
		}
		PerCohortValueReporter::UniquePtr pcv( new PerCohortValueReporter( names ) ) ;
		if( options.check( "-info-columns" )) {
			BOOST_FOREACH( std::string const& variable, options.get_values( "-info-columns" )) {
				pcv->add_variable( variable ) ;
			}
		}
		result.reset( pcv.release() ) ;
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

		{
			options.declare_group( "File handling options" ) ;
			options[ "-data" ]
				.set_description( "Specify the path of a file containing SNPTEST or MMM results to load." )
				.set_is_required()
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 1 )
				.set_maximum_multiplicity( 100 )
			;
			
			options[ "-info-columns" ]
				.set_description( "Specify extra columns in input files whose values will be considered as variables to be reported in the output."
				 	" Currently these must be columns of numerical data.  (A single wildcard character * at the start or end of the column name may"
					" be used to match any initial or terminal sequence of characters; but take care to escape this from the shell.)" )
					.set_takes_values_until_next_option()
					.set_minimum_multiplicity( 0 )
					.set_maximum_multiplicity( 100 )
				;
				
			options[ "-snp-match-fields" ]
				.set_description( "Use this option to specify a comma-separated list of SNP-identifying fields that should be used "
					"to match SNPs between cohorts.  Possible fields "
					"are \"position\", \"alleles\", \"rsid\", or \"snpid\"; you must always specify \"position\" as the first entry." )
				.set_takes_single_value()
				.set_default_value( "position,alleles" ) ;
		}
		{
			options.declare_group( "SNP exclusion options" ) ;
			options[ "-excl-snpids" ]
				.set_description( "Exclude all SNPs whose SNPID is in the given file(s) from the analysis.")
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-excl-rsids" ]
				.set_description( "Exclude all SNPs whose RSID is in the given file(s) from the analysis.")
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-snpids" ]
				.set_description( "Exclude all SNPs whose SNPID is not in the given file(s) from the analysis.")
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-rsids" ]
				.set_description( "Exclude all SNPs whose RSID is not in the given file(s) from the analysis.")
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-excl-positions" ]
				.set_description( "Exclude all SNPs whose position is in the given file(s) from the analysis. "
					"Positions should be in the form [chromosome]:[position] and separated by whitespace." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-positions" ]
				.set_description( "Exclude all SNPs whose position is not in the given file(s) from the analysis. "
					"Positions should be in the form [chromosome]:[position] and separated by whitespace." )
				.set_takes_values_until_next_option() 
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

			options[ "-excl-snpids-per-cohort" ]
				.set_description( "Exclude all SNPs whose SNPID is in the given file(s) from the corresponding cohort. "
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-excl-rsids-per-cohort" ]
				.set_description( "Exclude all SNPs whose RSID is in the given file(s) from the corresponding cohort. "
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-snpids-per-cohort" ]
				.set_description( "Exclude all SNPs whose SNPID is not in the given file(s) from the corresponding cohort."
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
			options[ "-incl-rsids-per-cohort" ]
				.set_description( "Exclude all SNPs whose RSID is not in the given file(s) from the corresponding cohort."
					"The i-th value should be a comma-separated list of files containing rsids to exclude from the i-th cohort." )
				.set_takes_values_until_next_option() 
				.set_maximum_multiplicity( 100 ) ;
		
			options[ "-min-info" ]
				.set_description( "Treat SNPs with info less than the given threshhold as missing." )
				.set_takes_values( 1 ) ;

			options[ "-min-maf" ]
				.set_description( "Treat SNPs with maf (in controls) less than the given threshhold as missing." )
				.set_takes_values( 1 ) ;
				
				
		}

		{
			options.declare_group( "Options for adjusting variants" ) ;
			options[ "-flip-alleles" ]
				.set_description( "Specify that, where alleles do not match between cohorts, bingwa will try to match them by flipping." )
			;
		}

		{
			options.declare_group( "Options affecting output" ) ;

			options[ "-o" ]
				.set_description( "Specify the path to the output file." )
				.set_is_required()
				.set_takes_single_value() ;

			options[ "-flat-file" ]
				.set_description( "Specify the output file should be a flat file, not a db." ) ;
			options[ "-flat-table" ]
				.set_description( "Output all results for this analysis to one table with variables in columns and variants in rows. "
					"This overrides the default db output style, which is in a normalised form with different variables on different rows." )
			;
			options[ "-analysis-name" ]
				.set_description( "Specify a name to label results from this analysis with" )
				.set_takes_single_value()
				.set_default_value( "bingwa analysis, started " + appcontext::get_current_time_as_string() ) ;
			options[ "-cohort-names" ]
				.set_description( "Specify a name to label results from this analysis with" )
				.set_takes_values_until_next_option() ;
			options[ "-table-name" ]
				.set_description( "Specify a name for the table to use when using -flat-table." )
				.set_takes_single_value() ;

			options[ "-log" ]
				.set_description( "Specify the path of a log file; all screen output will be copied to the file." )
				.set_takes_single_value() ;
				
			options.option_implies_option( "-table-name", "-flat-table" ) ;
		}

		{
			options.declare_group( "Analysis options" ) ;
		
			options[ "-no-meta-analysis" ]
				.set_description( "Don't do a fixed effect meta-analysis.  Instead, just match up SNPs and store per-cohort values." ) ;
				
			options[ "-simple-prior" ]
				.set_description( "Specify the between-cohort prior correlation, denoted r, giving a correlation matrix of the form\n"
					"   [ 1 r r  .. ]\n"
					"   [ r 1 r  .. ]\n"
					"   [ r r 1  .. ]\n"
					"   [       .   ]\n"
					"   [         . ]\n"
				)
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 ) ;

			options[ "-complex-prior" ]
				.set_description( "Specify the upper triangle of the prior correlation matrix for bayesian analysis.  For example, "
					"the value \"1,0.5,1\" specifies the matrix\n"
					"   [ 1    0.5 ]\n"
					"   [ 0.5  1   ]."
				)
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 ) ;

			options[ "-complex-prior-name" ]
				.set_description( "Specify the name of the complex models.  This option takes the same number of values as "
					"are given to -complex-prior option."
				)
				.set_takes_values_until_next_option()
				.set_minimum_multiplicity( 0 )
				.set_maximum_multiplicity( 1 )
			;

			options[ "-prior-sd" ]
				.set_description( "Specify the prior standard deviation for bayesian analysis." )
				.set_takes_values_until_next_option() ;

			options.option_excludes_option( "-no-meta-analysis", "-simple-prior" ) ;
			options.option_excludes_option( "-no-meta-analysis", "-complex-prior" ) ;
			options.option_implies_option( "-simple-prior", "-prior-sd" ) ;
			options.option_implies_option( "-complex-prior", "-prior-sd" ) ;
			options.option_implies_option( "-complex-prior", "-complex-prior-name" ) ;
		}
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
		m_snps( compare_fields ),
		m_flip_alleles_if_necessary( false )
	{
	}
	
	void set_flip_alleles( void ) {
		m_flip_alleles_if_necessary = true ;
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
				double frequency ;
				m_cohorts[i].get_info( snp_i, &info ) ;
				m_cohorts[i].get_frequency( snp_i, &frequency ) ;
				ui_context.logger() << "     " << m_cohorts[i].get_SNP( snp_i ) << " (frequency = " << frequency << ", info = " << info << ")\n";
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
			genfile::SNPIdentifyingData2 const& snp,
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
	struct SnpMatch {
	public:
		SnpMatch( std::size_t index_, bool flip_ ): index( index_ ), flip( flip_ ) {}
		SnpMatch(): index(0), flip( false ) {}
		SnpMatch( SnpMatch const& other ): index( other.index ), flip( other.flip ) {}
		SnpMatch& operator=( SnpMatch const& other ) { index = other.index ; flip = other.flip ; return *this ; }
	public:
		std::size_t index ;
		bool flip ;
	} ;
	typedef boost::optional< SnpMatch > OptionalSnpMatch ;
	typedef std::map<
		genfile::SNPIdentifyingData2,
		std::vector< OptionalSnpMatch >,
		genfile::SNPIdentifyingData2::CompareFields
	> SnpMap ;
	SnpMap m_snps ;
	typedef std::map< std::vector< bool >, std::size_t > CategoryCounts ;
	CategoryCounts m_category_counts ;
	bool m_flip_alleles_if_necessary ;
	
	ResultSignal m_result_signal ;

	struct DataGetter: public AmetComputation::DataGetter {
		DataGetter(
			boost::ptr_vector< FrequentistGenomeWideAssociationResults > const& cohorts,
			std::vector< OptionalSnpMatch >const& indices
		):
			m_cohorts( cohorts ),
			m_indices( indices )
		{}
		
		std::size_t get_number_of_cohorts() const { return m_cohorts.size() ; }
		
		void get_counts( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_counts( m_indices[i]->index, result ) ;
				if( m_indices[i]->flip ) {
					result->reverseInPlace() ;
				}
			}
		}
		void get_betas( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_betas( m_indices[i]->index, result ) ;
				if( m_indices[i]->flip ) {
					(*result) *= -1 ;
				}
			}
		}
		void get_ses( std::size_t i, Eigen::VectorXd* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_ses( m_indices[i]->index, result ) ;
			}
		}
		void get_pvalue( std::size_t i, double* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_pvalue( m_indices[i]->index, result ) ;
			}
		}
		void get_info( std::size_t i, double* result ) const {
			if( is_non_missing( i ) ) {
				m_cohorts[i].get_info( m_indices[i]->index, result ) ;
			}
		}
		bool is_non_missing( std::size_t i ) const {
			return( m_indices[i] ) ;
		}
		void get_variable( std::string const& variable, std::size_t i, double* result ) const {
			m_cohorts[i].get_variable( m_indices[i]->index, variable, result ) ;
		}

		private:
			boost::ptr_vector< FrequentistGenomeWideAssociationResults > const& m_cohorts ;
			std::vector< OptionalSnpMatch > const& m_indices ;
	} ;
private:
	void unsafe_setup( appcontext::UIContext& ui_context ) {
		link_data( ui_context ) ;
		categorise_by_missingness() ;
	}
	
	void add_SNP_callback( std::size_t cohort_i, std::size_t snp_i, genfile::SNPIdentifyingData2 const& snp ) {
		// Find the SNP that matches the given one (if it exists)
		std::pair< SnpMap::iterator, SnpMap::iterator > range = m_snps.equal_range( snp ) ;
		SnpMatch snp_match( snp_i, false ) ;

		if( range.second == range.first && m_flip_alleles_if_necessary ) {
			genfile::SNPIdentifyingData2 swapped_snp = snp ;
			swapped_snp.swap_alleles() ;
			range = m_snps.equal_range( swapped_snp ) ;
			snp_match.flip = ( range.second != range.first ) ;
		}
		
		if( range.second == range.first ) {
			// no match, so add this SNP.
			std::vector< OptionalSnpMatch >& snp_matches = m_snps[ snp ] ;
			snp_matches.resize( m_cohorts.size() ) ;
			snp_matches[ cohort_i ] = snp_match ;
		}
		else {
			// There is a match.  In case rsids differ, we combine them separated by commas.
			// First save the currently-stored data.
			genfile::SNPIdentifyingData2 stored_snp = range.first->first ;
			std::vector< OptionalSnpMatch > snp_matches = range.first->second ;

			m_snps.erase( range.first ) ;
			
			using genfile::string_utils::slice ;
			std::string const rsid_string = stored_snp.get_rsid() ;
			std::vector< slice > split_rsids = slice( rsid_string ).split( "," ) ;
			std::sort( split_rsids.begin(), split_rsids.end() ) ;
			if( !std::binary_search( split_rsids.begin(), split_rsids.end(), snp.get_rsid() ) ) {
				split_rsids.push_back( snp.get_rsid() ) ;
				std::sort( split_rsids.begin(), split_rsids.end() ) ;
				stored_snp.set_rsid( genfile::string_utils::join( split_rsids, "," ) ) ;
			}
			
			stored_snp.add_identifier( snp.get_rsid() ) ;
			snp.get_identifiers( boost::bind( &genfile::SNPIdentifyingData2::add_identifier, &stored_snp, _1 ) ) ;

			snp_matches[ cohort_i ] = snp_match ;
			m_snps[ stored_snp ] = snp_matches ;
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
	
	void categorise_by_missingness() {
		SnpMap::const_iterator snp_i = m_snps.begin(), end_i = m_snps.end() ;
		for( ; snp_i != end_i; ++snp_i ) {
			std::vector< OptionalSnpMatch > const& indices = snp_i->second ;
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
			appcontext::UIContext::ProgressContext progress_context = ui_context.get_progress_context( "Storing meta-analysis results" ) ;
			SnpMap::const_iterator snp_i = m_snps.begin() ;
			SnpMap::const_iterator const end_i = m_snps.end() ;
			for( std::size_t snp_index = 0; snp_i != end_i; ++snp_i, ++snp_index ) {
				std::vector< OptionalSnpMatch > const& indices = snp_i->second ;
				DataGetter data_getter( m_cohorts, indices ) ;
				for( std::size_t i = 0; i < m_computations.size(); ++i ) {
					m_computations[i](
						snp_i->first,
						data_getter,
						boost::bind(
							boost::ref( m_result_signal ),
							snp_i->first,
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
	{
		if( options().check( "-flip-alleles" )) {
			m_processor->set_flip_alleles() ;
		}
	}
	
	void run() {
		try {
			unsafe_run() ;
		}
		catch( genfile::InputError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::FileNotFoundError const& e ) {
			get_ui_context().logger() << "\nError: No file matching \"" << e.filespec() << "\" could be found.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}

	void unsafe_run() {
		load_data() ;

		m_processor->setup( get_ui_context() ) ;
		m_processor->summarise( get_ui_context() ) ;
		
		qcdb::Storage::SharedPtr storage ;
		if( options().check( "-flat-file" )) {
			storage = qcdb::FlatFileOutputter::create_shared(
				options().get< std::string >( "-o" ),
				options().get< std::string >( "-analysis-name" ),
				options().get_values_as_map()
			) ;
		}
		else if( options().check( "-flat-table" )) {
			qcdb::FlatTableDBOutputter::SharedPtr table_storage = qcdb::FlatTableDBOutputter::create_shared(
				options().get< std::string >( "-o" ),
				options().get< std::string >( "-analysis-name" ),
				options().get_values_as_map()
			) ;

			if( options().check( "-table-name" ) ) {
				table_storage->set_table_name( options().get< std::string >( "-table-name" )) ;
			}
			storage = table_storage ;
		}
		else {
			storage = snp_summary_component::DBOutputter::create_shared(
				options().get< std::string >( "-o" ),
				options().get< std::string >( "-analysis-name" ),
				options().get_values_as_map()
			) ;
		}

		m_processor->send_results_to(
			boost::bind(
				&qcdb::Storage::store_per_variant_data,
				storage,
				_1, _2, _3
			)
		) ;

		if( options().check( "-data" ) && options().get_values< std::string >( "-data" ).size() > 0 ) {
			m_processor->add_computation(
				"PerCohortValueReporter",
				AmetComputation::create( "PerCohortValueReporter", options() )
			) ;
			if( !options().check( "-no-meta-analysis" )) {
				m_processor->add_computation(
					"FixedEffectFrequentistMetaAnalysis",
					AmetComputation::create( "FixedEffectFrequentistMetaAnalysis", options() )
				) ;
				if( options().check( "-simple-prior" ) || options().check( "-complex-prior" )) {
					std::map< std::string, Eigen::MatrixXd > const priors = get_priors( options() ) ;
					assert( priors.size() > 0 ) ;
					std::map< std::string, Eigen::MatrixXd >::const_iterator i = priors.begin() ;
					std::map< std::string, Eigen::MatrixXd >::const_iterator const end_i = priors.end() ;
					for( ; i != end_i; ++i ) {
						ApproximateBayesianMetaAnalysis::UniquePtr computation(
							new ApproximateBayesianMetaAnalysis(
								i->first,
								i->second
							)
						) ;
						m_processor->add_computation(
							"ApproximateBayesianMetaAnalysis",
							AmetComputation::UniquePtr( computation.release() )
						) ;
					}
				}
			}
		}
		m_processor->process( get_ui_context() ) ;
		
		get_ui_context().logger() << "Finalising storage...\n" ;
		storage->finalise() ;
	}
	
	std::map< std::string, Eigen::MatrixXd > get_priors( appcontext::OptionProcessor const& options ) {
		std::map< std::string, Eigen::MatrixXd > result ;
		using genfile::string_utils::to_string ;
		using genfile::string_utils::to_repr ;
		using genfile::string_utils::split_and_strip_discarding_empty_entries ;
		int const N = options.get_values< std::string >( "-data" ).size() ;
		std::vector< std::string > sds = options.get_values< std::string >( "-prior-sd" ) ;
		
		std::size_t number_of_sds_used = 0 ;
		
		if( options.check( "-simple-prior" ) ) {
			std::vector< std::string > const rhos = options.get_values< std::string >( "-simple-prior" ) ;
			if( sds.size() < rhos.size() ) {
				throw genfile::BadArgumentError(
					"AmetProcessor::get_priors()",
					"-simple-prior \"" + options.get< std::string >( "-simple-prior" ) + "\""
				) ;
			}
			for( std::size_t i = 0; i < rhos.size(); ++i ) {
				double const rho = to_repr< double >( rhos[i] ) ;
				double const sd = to_repr< double >( sds[i] ) ;
				result[ "rho=" + rhos[i] + "/sd=" + sds[i] ] = get_prior_matrix( N, rho, sd ) ;
			}
			
			number_of_sds_used += rhos.size() ;
		}

		if( options.check( "-complex-prior" ) ) {
			std::vector< std::string > const matrix_specs = options.get_values< std::string >( "-complex-prior" ) ;
			std::vector< std::string > const model_names = options.get_values< std::string >( "-complex-prior-name" ) ;
			if( ( sds.size() - number_of_sds_used ) < matrix_specs.size() ) {
				throw genfile::BadArgumentError(
					"AmetProcessor::get_priors()",
					"-complex-prior \"" + options.get< std::string >( "-matrix-specs" ) + "\""
				) ;
			}
			if( model_names.size() != matrix_specs.size() ) {
				throw genfile::BadArgumentError(
					"AmetProcessor::get_priors()",
					"-complex-prior-name \"" + options.get< std::string >( "-complex-prior-name" ) + "\""
				) ;
			}
			for( std::size_t i = 0; i < matrix_specs.size(); ++i ) {
				double const sd = to_repr< double >( sds[i + number_of_sds_used] ) ;
				result[ model_names[i] + "/sd=" + sds[i + number_of_sds_used] ] = get_prior_matrix( 
					N,
					matrix_specs[i],
					sd
				) ;
			}
			number_of_sds_used += matrix_specs.size() ;
		}
		
		if( sds.size() != number_of_sds_used ) {
			throw genfile::BadArgumentError(
				"AmetProcessor::get_priors()",
				"-prior-sd \"" + options.get< std::string >( "-prior-sd" ) + "\""
			) ;
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
		std::vector< std::string > values = genfile::string_utils::split_and_strip( matrix_spec, "," ) ;
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
		boost::function< void ( genfile::SNPIdentifyingData2 const&, std::string const&, genfile::VariantEntry const& ) > 
		ResultCallback ;

	void load_data() {
		using genfile::string_utils::to_string ;

		std::vector< std::string > cohort_files = options().get_values< std::string >( "-data" ) ;
		for( std::size_t cohort_i = 0; cohort_i < cohort_files.size(); ++cohort_i ) {
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading scan results \"" + cohort_files[cohort_i] + "\"" ) ;

			SNPTESTResults::SNPResultCallback snp_callback ;
			if( cohort_files.size() == 1 ) {
				qcdb::Storage::SharedPtr raw_storage = snp_summary_component::DBOutputter::create_shared(
					options().get< std::string >( "-o" ),
					globals::program_name + " cohort " + to_string( cohort_i+1 ) + ": SNPTEST analysis: \"" + cohort_files[cohort_i] + "\"",
					snp_summary_component::DBOutputter::Metadata()
				) ;

				snp_callback = boost::bind(
					&qcdb::Storage::store_per_variant_data,
					raw_storage,
					_1, _2, _3
				) ;
			}

			SNPExclusionTestConjunction::UniquePtr test( new SNPExclusionTestConjunction() ) ;

			if( options().check_if_option_was_supplied_in_group( "SNP exclusion options" ) ) {
				genfile::CommonSNPFilter::UniquePtr snp_filter = get_snp_exclusion_filter( cohort_i ) ;
				if( snp_filter.get() ) {
					test->add_subtest( genfile::SNPIdentifyingDataTest::UniquePtr( snp_filter.release() ) ) ;
				}
				if( options().check( "-min-info" ) || options().check( "-min-maf" ) ) {
					ThreshholdingSummaryStatisticTest::UniquePtr subtest( new ThreshholdingSummaryStatisticTest() ) ;
					if( options().check( "-min-info" ) ) {
						subtest->set_inclusion_bounds( "info", options().get< double >( "-min-info" ), 2 ) ; // info should be no larger than 1, but I guess there may be numerical error.
					}
					if( options().check( "-min-maf" ) ) {
						subtest->set_inclusion_bounds( "maf", options().get< double >( "-min-maf" ), 1 ) ;
					}
					test->add_subtest( SNPExclusionTest::UniquePtr( subtest.release() ) ) ;
				}
			}

			FrequentistGenomeWideAssociationResults::UniquePtr results
				= FrequentistGenomeWideAssociationResults::create(
					genfile::wildcard::find_files_by_chromosome( cohort_files[cohort_i] ),
					options().check( "-info-columns" ) ? options().get_values< std::string >( "-info-columns" ) : std::vector< std::string >(),
					SNPExclusionTest::UniquePtr( test.release() ),
					snp_callback,
					progress_context
			) ;
			
			// summarise right now so as to see memory used.
			get_ui_context().logger() << "Cohort " << (cohort_i+1) << " summary: " << results->get_summary() << ".\n" ;
			
			m_processor->add_cohort( "cohort_" + to_string( cohort_i+1 ), results ) ;
		}
	}
	
	genfile::CommonSNPFilter::UniquePtr get_snp_exclusion_filter( std::size_t cohort_i ) const {
		genfile::CommonSNPFilter::UniquePtr snp_filter ;
		using genfile::string_utils::join ;
		using genfile::string_utils::split_and_strip_discarding_empty_entries ;

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
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-snpids" ) ;
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

			if( options().check_if_option_was_supplied( "-excl-positions" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-positions" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::Positions
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-positions" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-positions" ) ;
				BOOST_FOREACH( std::string const& filename, files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::Positions
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-snpids-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-snpids-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "AmetApplication::get_snp_exclusion_filter()", "-excl-snpids-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::SNPIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-snpids-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-snpids-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "AmetApplication::get_snp_exclusion_filter()", "-incl-snpids-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::SNPIDs
					) ;
				}
			}


			if( options().check_if_option_was_supplied( "-excl-rsids-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-excl-rsids-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "AmetApplication::get_snp_exclusion_filter()", "-excl-rsids-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-rsids-per-cohort" )) {
				std::vector< std::string > files = options().get_values< std::string > ( "-incl-rsids-per-cohort" ) ;
				if( files.size() != options().get_values< std::string >( "-data" ).size() ) {
					throw genfile::BadArgumentError( "AmetApplication::get_snp_exclusion_filter()", "-incl-rsids-per-cohort=\"" + join( files, " " ) + "\"" ) ;
				}
				std::vector< std::string > this_cohort_files = split_and_strip_discarding_empty_entries( files[ cohort_i ], ",", " \t" ) ;
				BOOST_FOREACH( std::string const& filename, this_cohort_files ) {
					snp_filter->include_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-snps-matching" )) {
				std::string const it = options().get< std::string > ( "-excl-snps-matching" ) ;
				std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( it, ",", " \t" ) ;
				BOOST_FOREACH( std::string const& spec, specs ) {
					snp_filter->exclude_snps_matching(
						spec
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-incl-snps-matching" )) {
				std::string const it = options().get< std::string > ( "-incl-snps-matching" ) ;
				std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( it, ",", " \t" ) ;
				BOOST_FOREACH( std::string const& spec, specs ) {
					snp_filter->include_snps_matching(
						spec
					) ;
				}
			}
			
			if( options().check_if_option_was_supplied( "-incl-range" )) {
				std::vector< std::string > const specs = options().get_values< std::string >( "-incl-range" ) ;
				for ( std::size_t i = 0; i < specs.size(); ++i ) {
					snp_filter->include_snps_in_range(
						genfile::GenomePositionRange::parse( specs[i] )
					) ;
				}
			}

			if( options().check_if_option_was_supplied( "-excl-range" )) {
				std::vector< std::string > specs = options().get_values< std::string >( "-excl-range" ) ;
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
		AmetApplication app( argc, argv ) ;	
		app.run() ;
		app.post_summarise() ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
