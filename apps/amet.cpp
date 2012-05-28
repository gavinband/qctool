
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <sstream>
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
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/CommonSNPFilter.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"

namespace globals {
	std::string const program_name = "amet" ;
	std::string const program_version = "0.1" ;
}

struct FrequentistGenomeWideAssociationResults: public boost::noncopyable {
	typedef std::auto_ptr< FrequentistGenomeWideAssociationResults > UniquePtr ;
	virtual ~FrequentistGenomeWideAssociationResults() {}
	typedef boost::function< void ( std::size_t i, genfile::SNPIdentifyingData const& snp ) > SNPCallback ;
	virtual void get_SNPs( SNPCallback ) const = 0 ;
	virtual void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ;
	virtual void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ; 
	virtual void get_pvalue( std::size_t snp_i, double* result ) const = 0 ;
	virtual void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const = 0 ;
	virtual std::string get_summary( std::string const& prefix = "", std::size_t target_column = 80 ) const = 0;
} ;

struct SNPTESTResults: public FrequentistGenomeWideAssociationResults {
	typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;
	typedef
		boost::function< void ( std::size_t, genfile::SNPIdentifyingData const&, std::string const&, std::string const&, genfile::VariantEntry const& ) > 
		SNPResultCallback ;
	
	SNPTESTResults( std::vector< genfile::wildcard::FilenameMatch > const& filenames, SNPResultCallback callback = SNPResultCallback(), ProgressCallback progress_callback = ProgressCallback() )
	{
		setup( filenames, callback, progress_callback ) ;
	}

	void get_SNPs( SNPCallback callback ) const {
		for( std::size_t i = 0; i < m_snps.size(); ++i ) {
			callback( i, m_snps[i] ) ;
		}
	}

	void get_betas( std::size_t snp_i, Eigen::VectorXd* result ) const {
		*result = m_betas.row( snp_i ) ;
	}
	void get_ses( std::size_t snp_i, Eigen::VectorXd* result ) const {
		*result = m_ses.row( snp_i ) ;
	}
	void get_pvalue( std::size_t snp_i, double* result ) const {
		*result = m_pvalues( snp_i ) ;
	}
	void get_counts( std::size_t snp_i, Eigen::VectorXd* result ) const {
		*result = m_sample_counts.row( snp_i ) ;
	}

	std::string get_summary( std::string const& prefix, std::size_t target_column ) const {
		using genfile::string_utils::to_string ;
		std::string result = prefix + "SNPTESTResults object holding results for " + to_string( m_snps.size() ) + " SNPs.  The first few are:\n" ;
		for( std::size_t i = 0; i < std::min( std::size_t( 5 ), m_snps.size() ); ++i ) {
			result += to_string( m_snps[i] ) + ": beta=" + to_string( m_betas.row(i) ) + ", se=" + to_string( m_ses.row(i) ) + ", pvalue=" + to_string( m_pvalues(i) ) + "\n" ;
		}
		return result ;
	}

private:
	std::vector< genfile::SNPIdentifyingData > m_snps ;
	Eigen::MatrixXd m_betas ;
	Eigen::MatrixXd m_ses ;
	Eigen::VectorXd m_pvalues ;
	Eigen::VectorXd m_info ;
	Eigen::MatrixXd m_sample_counts ;
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
	
	void setup( std::size_t const filename_i, std::size_t const number_of_files, genfile::wildcard::FilenameMatch const& filename, SNPResultCallback callback,ProgressCallback progress_callback ) {
		statfile::BuiltInTypeStatSource::UniquePtr source( statfile::BuiltInTypeStatSource::open( filename.filename() )) ;

		ColumnMap column_map = get_columns_to_store( *source ) ;
		
		int degrees_of_freedom = ( column_map.left.find( "_beta_2" ) == column_map.left.end() ) ? 1 : 2 ;
		
		m_snps.resize( source->number_of_rows() ) ;
		m_betas.resize( source->number_of_rows(), degrees_of_freedom ) ;
		m_ses.resize( source->number_of_rows(), degrees_of_freedom ) ;
		m_pvalues.resize( source->number_of_rows() ) ;
		m_info.resize( source->number_of_rows() ) ;
		m_sample_counts.resize( source->number_of_rows(), 4 ) ;

		genfile::SNPIdentifyingData snp ;
		for(
			std::size_t snp_index = 0 ;
			(*source) >> snp.SNPID() >> snp.rsid() >> snp.position().chromosome() >> snp.position().position() >> snp.first_allele() >> snp.second_allele();
			++snp_index, (*source) >> statfile::ignore_all()
		) {
			m_snps[ snp_index ] = snp ;
	
			// Just read the values we need
			ColumnMap::right_const_iterator i = column_map.right.begin(), end_i = column_map.right.end() ;
			double value ;
			
			for( ; i != end_i; ++i ) {
				if( 0 ) { //callback ) {
					// We have to read all the columns, even though we'll ignore them.
					while( source->current_column() < source->number_of_columns() && ( i == end_i || source->current_column() < i->first ) ) {
						(*source) >> value ;
						callback( snp_index, snp, "SNPTESTResults", source->name_of_column( source->current_column() - 1 ), value ) ;
					}
				}
				else {
					(*source) >> statfile::ignore( i->first - source->current_column() ) ;
				}
				(*source) >> value ;
				if( callback ) {
					callback( snp_index, snp, "SNPTESTResults", source->name_of_column( source->current_column() - 1 ), value ) ;
				}
				store_value( snp_index, i->second, value ) ;
			}
			
			// read any remaining columns if necessary.
			if( callback ) {
				while( source->current_column() < source->number_of_columns() && ( i == end_i || source->current_column() < i->first ) ) {
					(*source) >> value ;
					callback( snp_index, snp, "SNPTESTResults", source->name_of_column( source->current_column() - 1 ), value ) ;
				}
			}
			
			if( progress_callback ) {
				progress_callback( 100 * ( filename_i + double( source->number_of_rows_read() ) / source->number_of_rows() ), 100 * number_of_files ) ;
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

struct FilteringGenomeWideAssociationResults: public FrequentistGenomeWideAssociationResults {
	FilteringGenomeWideAssociationResults( FrequentistGenomeWideAssociationResults::UniquePtr source, genfile::SNPIdentifyingDataTest::UniquePtr test ):
		m_source( source ),
		m_test( test )
	{
		assert( m_test.get() ) ;
		link_indices() ;
	}

	// typedef boost::function< void ( std::size_t i, genfile::SNPIdentifyingData const& snp ) > SNPCallback ;

	void get_SNPs( SNPCallback callback ) const {
		get_filtered_SNPs( callback ) ;
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
	std::string get_summary( std::string const& prefix = "", std::size_t target_column = 80 ) const {
		return prefix + "FilteringGenomeWideAssociationResults object holding " + genfile::string_utils::to_string( m_links.size() ) + " SNPs, from: " + m_source->get_summary( "", target_column ) ;
	}
	
private:
	FrequentistGenomeWideAssociationResults::UniquePtr m_source ;
	genfile::SNPIdentifyingDataTest::UniquePtr m_test ;
	std::vector< std::size_t > m_links ;
private:
	
	struct SNPFilterCallback {
		SNPFilterCallback( genfile::SNPIdentifyingDataTest const& test, SNPCallback callback = SNPCallback() ):
			m_test( test ),
			m_callback( callback ),
			m_i( 0 )
		{
		} ;
		
		void operator()( std::size_t i, genfile::SNPIdentifyingData const& snp ) {
			if( m_test( snp ) ) {
				m_links.resize( m_i + 1 ) ;
				m_links[ m_i ] = i ;
				if( m_callback ) {
					m_callback( m_i, snp ) ;
				}
				m_i++ ;
			}
			else {
			}
		}

		std::vector< std::size_t > const& get_links() const { return m_links ; }

	private:
		genfile::SNPIdentifyingDataTest const& m_test ;
		SNPCallback m_callback ;
		std::size_t m_i ;
		std::vector< std::size_t > m_links ;
	} ;
	
	// populate the link between our filtered SNP indices and the indices in the base source.
	void link_indices() {
		assert( m_test.get() ) ;
		SNPFilterCallback filtering_callback( *m_test ) ;
		m_source->get_SNPs( boost::ref( filtering_callback ) ) ;
		m_links = filtering_callback.get_links() ;
	}
	
	void get_filtered_SNPs( SNPCallback callback ) const {
		assert( m_test.get() ) ;
		SNPFilterCallback filtering_callback( *m_test, callback ) ;
		m_source->get_SNPs( filtering_callback ) ;
	}
} ;

struct AmetComputation: public boost::noncopyable {
	typedef std::auto_ptr< AmetComputation > UniquePtr ;
	static UniquePtr create( std::string const& name, appcontext::OptionProcessor const& ) ;
	virtual ~AmetComputation() {}
	typedef genfile::SNPIdentifyingData SNPIdentifyingData ;
	typedef boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;

	struct DataGetter: public boost::noncopyable {
		virtual ~DataGetter() {} ;
		virtual std::size_t get_number_of_cohorts() const = 0 ;
		virtual bool is_non_missing( std::size_t i ) const = 0 ;
		virtual void get_counts( std::size_t, Eigen::VectorXd* result ) const = 0 ;
		virtual void get_betas( std::size_t i, Eigen::VectorXd* result ) const = 0 ;
		virtual void get_ses( std::size_t i, Eigen::VectorXd* result  ) const = 0 ;
		virtual void get_pvalue( std::size_t i, double* result ) const = 0 ;
	} ;

	virtual void operator()(
		SNPIdentifyingData const&,
		DataGetter const& data_getter,
		ResultCallback callback
	) = 0 ;
	virtual std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const = 0 ;
	virtual std::string get_spec() const = 0 ;
} ;

struct PerCohortValueReporter: public AmetComputation {
	void operator()(
		SNPIdentifyingData const&,
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
				data_getter.get_betas( i, &betas ) ;
				data_getter.get_ses( i, &ses ) ;
				data_getter.get_counts( i, &counts ) ;
				data_getter.get_pvalue( i, &pvalue ) ;
		
				assert( counts.size() == 4 ) ;
				using genfile::string_utils::to_string ;
				std::string prefix = "cohort_" + to_string( i + 1 ) + ":";
				callback( prefix + "AA", counts(0) ) ;
				callback( prefix + "AB", counts(0) ) ;
				callback( prefix + "BB", counts(0) ) ;
				callback( prefix + "NULL", counts(0) ) ;
				
				assert( betas.size() == ses.size() ) ;
				for( int j = 0; j < betas.size(); ++j ) {
					callback( prefix + "beta_" + to_string( j+1 ), betas(j) ) ;
					callback( prefix + "se_" + to_string( j+1 ), betas(j) ) ;
				}
				callback( prefix + "pvalue", pvalue ) ;
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
		SNPIdentifyingData const& snp,
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
			callback( "FixedEffectFrequentistMetaAnalysis/study_betas", impl::serialise( betas ) ) ;
			callback( "FixedEffectFrequentistMetaAnalysis/study_ses", impl::serialise( ses ) ) ;
			
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
		Eigen::MatrixXd const& sigma
	):
		m_sigma( sigma )
	{
		assert( m_sigma.rows() == m_sigma.cols() ) ;
	}

	void operator()(
		SNPIdentifyingData const& snp,
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
			callback( "ApproximateBayesianMetaAnalysis/bf", genfile::MissingValue() ) ;
			return ;
		}
		else {
			callback( "ApproximateBayesianMetaAnalysis/study_betas", impl::serialise( betas ) ) ;
			callback( "ApproximateBayesianMetaAnalysis/study_ses", impl::serialise( ses ) ) ; 
			
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
			
			// I hope LDLT copes with noninvertible matrices.
			// Maybe it doesn't...but let's find out.
			Eigen::LDLT< Eigen::MatrixXd > Vsolver( V ) ;
			Eigen::LDLT< Eigen::MatrixXd > V_plus_prior_solver( V + prior ) ;
			
			Eigen::VectorXd exponent = betas.transpose() * ( Vsolver.solve( betas ) + V_plus_prior_solver.solve( betas ) ) ;
			
			assert( exponent.size() == 1 ) ;
			
			double const constant = Vsolver.vectorD().prod() / V_plus_prior_solver.vectorD().prod() ;
			double const result = constant * std::exp( 0.5 * exponent(0) ) ;
			
			callback( "ApproximateBayesianMetaAnalysis/bf", result ) ;
		}
	}

	std::string get_spec() const {
		return "ApproximateBayesianMetaAnalysis with prior:\n" + genfile::string_utils::to_string( m_sigma ) ;
	}

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}
private:
	Eigen::MatrixXd const m_sigma ;
} ;


AmetComputation::UniquePtr AmetComputation::create( std::string const& name, appcontext::OptionProcessor const& options ) {
	AmetComputation::UniquePtr result ;
	if( name == "FixedEffectFrequentistMetaAnalysis" ) {
		result.reset( new FixedEffectFrequentistMetaAnalysis() ) ;
	}
	else if( name == "ApproximateBayesianMetaAnalysis" ) {
		std::vector< std::string > values = genfile::string_utils::split_and_strip( options.get< std::string >( "-prior-correlation" ), ",", " \t" ) ;
		
		double n = std::floor( std::sqrt( values.size() * 2 ) ) ;
		
		if( values.size() != (( n * (n+1) ) / 2 ) ) {
			throw genfile::BadArgumentError( "AmetComputation::create()", "-prior-correlation=\"" + options.get< std::string >( "-prior-correlation" ) + "\"" ) ;
		}
		
		Eigen::MatrixXd prior( n, n ) ;
		{
			int index = 0 ;
			for( int i = 0; i < n; ++i ) {
				for( int j = 0; j < n; ++j ) {
					prior( i, j )
						= ( j >= i ) ?
							genfile::string_utils::to_repr< double >( values[ index++ ] )
							:
							prior( j, i ) ;
				}
			}
		}
		
		prior *= options.get< double >( "-prior-variance" ) ;
		result.reset( new ApproximateBayesianMetaAnalysis( prior )) ;
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
		
		options[ "-excl-rsids" ]
			.set_description( "Specify a file containing a list of rsids to exclude from the analysis. "
				"If this option is specified it must occur the same number of times as the -snptest option." )
			.set_takes_values( 1 )
			.set_minimum_multiplicity( 0 )
			.set_maximum_multiplicity( 100 )
		;
		
		options[ "-prior-correlation" ]
			.set_description( "Specify the upper triangle of the prior correlation matrix for bayesian analysis." )
			.set_takes_single_value() ;

		options[ "-prior-variance" ]
			.set_description( "Specify the prior variance for bayesian analysis." )
			.set_takes_single_value() ;
		
		options.option_implies_option( "-prior-variance", "-prior-correlation" ) ;
		options.option_implies_option( "-prior-correlation", "-prior-variance" ) ;

		options[ "-omit-raw-results" ]
			.set_description( "Indicate that " + globals::program_name + " should not store raw results from the input files in the database." ) ;
		
		options[ "-o" ]
			.set_description( "Specify the path to a file in which results will be placed." )
			.set_is_required()
			.set_takes_single_value() ;
		
		options[ "-analysis-name" ]
			.set_description( "Specify a name to label results from this analysis with" )
			.set_takes_single_value()
			.set_default_value( "amet analysis, started " + appcontext::get_current_time_as_string() ) ;
		
		options[ "-snp-match-fields" ]
			.set_description( "Set the fields used to describe SNPs as equal." )
			.set_takes_single_value()
			.set_default_value( "position,alleles" ) ;

		options[ "-log" ]
			.set_description( "Specify the path of the log file." )
			.set_takes_single_value() ;
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
	static UniquePtr create( genfile::SNPIdentifyingData::CompareFields const& compare_fields ) {
		return UniquePtr( new AmetProcessor( compare_fields ) ) ;
	}
	
	AmetProcessor( genfile::SNPIdentifyingData::CompareFields const& compare_fields ):
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
			ui_context.logger() << m_cohorts[i].get_summary( " - " ) ;
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
				ui_context.logger() << std::setw(12) << i->first[j] ;
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

	typedef boost::signals2::signal< void ( std::size_t index, genfile::SNPIdentifyingData const& snp, std::string const& computation_name, std::string const& value_name, genfile::VariantEntry const& value ) > ResultSignal ;

	void send_results_to( ResultSignal::slot_type callback ) {
		m_result_signal.connect( callback ) ;
	}
	
private:
	std::vector< std::string > m_cohort_names ;
	boost::ptr_vector< FrequentistGenomeWideAssociationResults > m_cohorts ;
	boost::ptr_vector< AmetComputation > m_computations ;
	typedef std::map< genfile::SNPIdentifyingData, std::vector< boost::optional< std::size_t > >, genfile::SNPIdentifyingData::CompareFields > SnpMap ;
	SnpMap m_snps ;
	typedef std::map< std::vector< bool >, std::size_t > CategoryCounts ;
	CategoryCounts m_category_counts ;
	
	ResultSignal m_result_signal ;

	struct DataGetter: public AmetComputation::DataGetter {
		DataGetter( boost::ptr_vector< FrequentistGenomeWideAssociationResults > const& cohorts, std::vector< boost::optional< std::size_t > >const& indices ):
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
	
	void add_SNP_callback( std::size_t cohort_i, std::size_t snp_i, genfile::SNPIdentifyingData const& snp ) {
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
			genfile::SNPIdentifyingData stored_snp = range.first->first ;
			std::vector< boost::optional< std::size_t > > snp_indices = range.first->second ;
			m_snps.erase( range.first ) ;
			if( stored_snp.get_rsid() != snp.get_rsid() ) {
				stored_snp.rsid() += "," + snp.get_rsid() ;
			}
			if( stored_snp.get_SNPID() != snp.get_SNPID() ) {
				stored_snp.SNPID() += "," + snp.get_SNPID() ;
			}
			if( stored_snp.get_first_allele() != snp.get_first_allele() ) {
				stored_snp.first_allele() += "," + snp.get_first_allele() ;
			}
			if( stored_snp.get_second_allele() != snp.get_second_allele() ) {
				stored_snp.second_allele() += "," + snp.get_second_allele() ;
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
		
		appcontext::UIContext::ProgressContext progress_context = ui_context.get_progress_context( "Storing meta-analysis results" ) ;
		
		SnpMap::const_iterator snp_i = m_snps.begin(), end_i = m_snps.end() ;
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
		m_processor( AmetProcessor::create( genfile::SNPIdentifyingData::CompareFields( options().get< std::string > ( "-snp-match-fields" )) ) )
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
				"FixedEffectFrequentistMetaAnalysis",
				AmetComputation::create( "FixedEffectFrequentistMetaAnalysis", options() )
			) ;
			if( options().check( "-prior-correlation" )) {
				m_processor->add_computation(
					"ApproximateBayesianMetaAnalysis",
					AmetComputation::create( "ApproximateBayesianMetaAnalysis", options() )
				) ;
			}
		}
		m_processor->process( get_ui_context() ) ;
	}
	
	typedef
		boost::function< void ( std::size_t, genfile::SNPIdentifyingData const&, std::string const&, std::string const&, genfile::VariantEntry const& ) > 
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

			results.reset(
				new SNPTESTResults(
					genfile::wildcard::find_files_by_chromosome( cohort_files[i] ),
					snp_callback,
					progress_context
				)
			) ;
			
			if( options().check( "-excl-rsids" ) ) {
				genfile::CommonSNPFilter::UniquePtr snp_filter ( new genfile::CommonSNPFilter ) ;

				std::vector< std::string > exclusion_files = options().get_values< std::string > ( "-excl-rsids" ) ;
				assert( exclusion_files.size() == cohort_files.size() ) ;
				
				std::vector< std::string > cohort_exclusion_files = genfile::string_utils::split_and_strip_discarding_empty_entries( exclusion_files[i], ",", " \t" ) ;
				foreach( std::string const& filename, cohort_exclusion_files ) {
					snp_filter->exclude_snps_in_file(
						filename,
						genfile::CommonSNPFilter::RSIDs
					) ;
				}
			
				results.reset(
					new FilteringGenomeWideAssociationResults( results, genfile::SNPIdentifyingDataTest::UniquePtr( snp_filter.release() ) )
				) ;
			}
			
			m_processor->add_cohort( "cohort_" + to_string( i+1 ), results ) ;
		}
	}
	
private:
	AmetProcessor::UniquePtr m_processor ;
} ;

int main( int argc, char **argv ) {
	try {
		AmetApplication app( argc, argv ) ;	
		app.run() ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
