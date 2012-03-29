#include <string>
#include <memory>
#include <set>
#include <map>
#include <utility>
#include <boost/bind.hpp>
#include <boost/signals2/signal.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/regex.hpp>
#include <boost/math/distributions/normal.hpp>
#include <Eigen/Core>
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/ProgramFlow.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/utility.hpp"
#include "genfile/SNPIdentifyingData.hpp"
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
	SNPTESTResults( std::vector< genfile::wildcard::FilenameMatch > const& filenames, ProgressCallback progress_callback = ProgressCallback() )
	{
		setup( filenames, progress_callback ) ;
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
		
		void setup( std::vector< genfile::wildcard::FilenameMatch > const& filenames, ProgressCallback progress_callback = ProgressCallback() ) {
			progress_callback( 0, filenames.size() ) ;
			for( std::size_t i = 0; i < filenames.size(); ++i ) {
				setup( filenames[i] ) ;
				progress_callback( i+1, filenames.size() ) ;
			}
		}
		
		void setup( genfile::wildcard::FilenameMatch const& filename ) {
			statfile::BuiltInTypeStatSource::UniquePtr source( statfile::BuiltInTypeStatSource::open( filename.filename() )) ;
			std::map< std::string, std::size_t > columns_by_name ;
			std::map< std::size_t, std::string > columns_by_index ;
			columns_by_name[ "_beta_1" ] = 0 ;
			columns_by_name[ "_se_1" ] = 0 ;
			columns_by_name[ "_beta_2" ] = 0 ;
			columns_by_name[ "_se_2" ] = 0 ;
			columns_by_name[ "_pvalue" ] = 0 ;
			columns_by_name[ "info" ] = 0 ;
			columns_by_name[ "all_AA" ] = 0 ;
			columns_by_name[ "all_AB" ] = 0 ;
			columns_by_name[ "all_BB" ] = 0 ;
			columns_by_name[ "all_NULL" ] = 0 ;
			for( std::size_t i = 0; i < source->number_of_columns(); ++i ) {
				std::string name = source->name_of_column( i ) ;
				for( std::map< std::string, std::size_t >::iterator j = columns_by_name.begin(); j != columns_by_name.end(); ++j ) {
					if( name.size() >= j->first.size() && name.compare( name.size() - j->first.size(), j->first.size(), j->first ) == 0 ) {
						j->second = i ;
						columns_by_index[i] = j->first ;
						
						//std::cerr << "Column \"" << j->first << "\" found at position " << i << ".\n" ;
					}
				}
			}
		
			// Figure out if this is a general test.
			int degrees_of_freedom ;
			if( columns_by_name[ "_beta_2" ] > 0 ) {
				degrees_of_freedom = 2 ;
			} else {
				columns_by_name.erase( "_beta_2" ) ;
				columns_by_name.erase( "_se_2" ) ;
				degrees_of_freedom = 1; 
			}
		
			m_snps.resize( source->number_of_rows() ) ;
			m_betas.resize( source->number_of_rows(), degrees_of_freedom ) ;
			m_ses.resize( source->number_of_rows(), degrees_of_freedom ) ;
			m_pvalues.resize( source->number_of_rows() ) ;
			m_info.resize( source->number_of_rows() ) ;
			m_sample_counts.resize( source->number_of_rows(), 4 ) ;

			using genfile::string_utils::to_repr ;
			
			genfile::SNPIdentifyingData snp ;
			for(
				std::size_t index = 0 ;
				(*source) >> snp.SNPID() >> snp.rsid() >> snp.position().chromosome() >> snp.position().position() >> snp.first_allele() >> snp.second_allele();
				++index, (*source) >> statfile::ignore_all()
			) {
				m_snps[ index ] = snp ;
		
				std::map< std::size_t, std::string >::const_iterator
					i = columns_by_index.begin(),
					end_i = columns_by_index.end() ;
				for( ; i != end_i; ++i ) {
					std::string value ;
					(*source) >> statfile::ignore( i->first - source->current_column() ) >> value ;
					
					if( i->second == "_beta_1" ) {
						m_betas( index, 0 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "_se_1" ) {
						m_ses( index, 0 ) = to_repr< double >( value ) ;
					}
					if( i->second == "_beta_2" ) {
						m_betas( index, 1 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "_se_2" ) {
						m_ses( index, 1 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "_pvalue" ) {
						m_pvalues( index ) = to_repr< double >( value ) ;
					}
					else if( i->second == "info" ) {
						m_info( index ) = to_repr< double >( value ) ;
					}
					else if( i->second == "all_AA" ) {
						m_sample_counts( index, 0 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "all_AB" ) {
						m_sample_counts( index, 1 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "all_BB" ) {
						m_sample_counts( index, 2 ) = to_repr< double >( value ) ;
					}
					else if( i->second == "all_NULL" ) {
						m_sample_counts( index, 3 ) = to_repr< double >( value ) ;
					}
				}
			}
		}
} ;

struct AmetComputation: public boost::noncopyable {
	typedef std::auto_ptr< AmetComputation > UniquePtr ;
	static UniquePtr create( std::string const& name ) ;
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

struct FixedEffectFrequentistMetaAnalysis: public AmetComputation {
	void operator()(
		SNPIdentifyingData const&,
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
		
		get_data( data_getter, betas, ses, non_missingness ) ;
		
		Eigen::VectorXd inverse_variances = ( ses.array() * ses.array() ).inverse() ;
		for( int i = 0; i < N; ++i ) {
			if( non_missingness( i ) == 0.0 ) {
				inverse_variances( i ) = 0 ;
				betas( i ) = 0 ;
				ses( i ) = 0 ;
			}
		}
		double const meta_beta = ( inverse_variances.array() * betas.array() ).sum() / inverse_variances.sum() ;
		double const meta_se = 1.0 / inverse_variances.sum() ;

		callback( "fixed_effect_meta_beta", meta_beta ) ;
		callback( "fixed_effect_meta_se", meta_se ) ;

		if( meta_se > 0 && meta_se != std::numeric_limits< double >::infinity() ) {
			typedef boost::math::normal NormalDistribution ;
			NormalDistribution normal( 0, meta_se ) ;
			// P-value is the mass under both tails of the normal distribution larger than |meta_beta|
			double const pvalue = 2.0 * boost::math::cdf( boost::math::complement( normal, std::abs( meta_beta ) ) ) ;
			callback( "fixed_effect_meta_pvalue", pvalue ) ;
		}
	}

	std::string get_spec() const {
		return "FixedEffectFrequentistMetaAnalysis" ;
	}

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
		return prefix + get_spec() ;
	}
	
private:
	
	void get_data( DataGetter const& data_getter, Eigen::VectorXd& betas, Eigen::VectorXd& ses, Eigen::VectorXd& non_missingness ) {
		std::size_t N = betas.size() ;
		for( std::size_t i = 0; i < N; ++i ) {
			Eigen::VectorXd data ;
			data_getter.get_betas( i, &data ) ;
			assert( data.size() == 1 ) ;
			betas(i) = data(0) ;
			data_getter.get_ses( i, &data ) ;
			assert( data.size() == 1 ) ;
			ses(i) = data(0) ;

			non_missingness( i ) = data_getter.is_non_missing( i ) ? 1.0 : 0.0 ;

			if( betas(i) != betas(i) || ses(i) != ses(i) ) {
				non_missingness(i) = 0 ;
			}
		}
	}
} ;

AmetComputation::UniquePtr AmetComputation::create( std::string const& name ) {
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
		construct_computations() ;
	}
	
	void add_SNP_callback( std::size_t cohort_i, std::size_t snp_i, genfile::SNPIdentifyingData const& snp ) {
		std::vector< boost::optional< std::size_t > >& snp_indices = m_snps[ snp ] ;
		snp_indices.resize( m_cohorts.size() ) ;
		snp_indices[ cohort_i ] = snp_i ;
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

	void construct_computations() {
		m_computations.push_back( AmetComputation::create( "FixedEffectFrequentistMetaAnalysis" )) ;
		m_computations.push_back( AmetComputation::create( "PerCohortValueReporter" )) ;
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
		
		impl::DBOutputter::SharedPtr outputter = impl::DBOutputter::create_shared(
			options().get< std::string >( "-o" ),
			options().get< std::string >( "-analysis-name" ),
			options().get_values_as_map()
		) ;
		
		m_processor->send_results_to(
			boost::bind(
				&impl::DBOutputter::operator(),
				outputter,
				_1, _2, _3, _4, _5
			)
		) ;

		m_processor->process( get_ui_context() ) ;
	}
	
	void load_data() {
		using genfile::string_utils::to_string ;
		std::vector< std::string > cohort_files = options().get_values< std::string >( "-snptest" ) ;
		for( std::size_t i = 0; i < cohort_files.size(); ++i ) {
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading SNPTEST results \"" + cohort_files[i] + "\"" ) ;
			
			FrequentistGenomeWideAssociationResults::UniquePtr results( new SNPTESTResults( genfile::wildcard::find_files_by_chromosome( cohort_files[i] ), progress_context ) ) ;
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
