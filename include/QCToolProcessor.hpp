#ifndef QCTOOL_PROCESSOR_HPP
#define QCTOOL_PROCESSOR_HPP

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <memory>
#include <numeric>
#include <boost/bind.hpp>

#include "Timer.hpp"
#include "GenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "RowCondition.hpp"
#include "GenRowStatistics.hpp"
#include "SampleRowStatistics.hpp"
#include "ObjectSource.hpp"
#include "SimpleFileObjectSource.hpp"
#include "SimpleFileObjectSink.hpp"
#include "OstreamTee.hpp"
#include "QCToolContext.hpp"

struct QCToolProcessorException: public std::exception
{
	char const* what() const throw() {return "QCToolProcessorException" ; }
} ;

struct GenAndSampleFileMismatchException: public QCToolProcessorException
{
	GenAndSampleFileMismatchException( GenRowIdentifyingData const& data, std::size_t actual_number_of_samples, std::size_t expected_number_of_samples )
		: m_expected_number_of_samples( expected_number_of_samples ),
		m_actual_number_of_samples( actual_number_of_samples ),
		m_data( data )
	{}

	~GenAndSampleFileMismatchException() throw() {}
	
	char const* what() const throw() {return "GenAndSampleFileMismatchException" ; }

	std::size_t expected_number_of_samples() const { return m_expected_number_of_samples ; }
	std::size_t actual_number_of_samples() const { return m_actual_number_of_samples ; }
	GenRowIdentifyingData data() const { return m_data ; }
private:
	std::size_t m_expected_number_of_samples, m_actual_number_of_samples ;
	GenRowIdentifyingData m_data ;
} ;

struct QCToolProcessor
{
public:
	QCToolProcessor( QCToolContext & context )
		: m_context( context )
	{
	}

	void process() {
		try {
			unsafe_process() ;
		}
		catch( StatisticNotFoundException const& e ) {
			std::cerr << "!! ERROR: " << e << ".\n" ;
			std::cerr << "Note: required statistics must be added using -statistics.\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( GToolException const& e) {
			std::cerr << "!! ERROR: " << e << ".\n" ;
			throw HaltProgramWithReturnCode( -1 ) ;
		}
		catch( GenAndSampleFileMismatchException const e ) {
			std::cerr << "!! ERROR: " << e.what() << ":\n"
				<< "I think the SNP identified as "
				<< e.data()
				<< " has the wrong number of samples.\n"
				<< "  (" << e.actual_number_of_samples() << " instead of " << e.expected_number_of_samples() << ").\n" ; 
			throw HaltProgramWithReturnCode( -1 ) ;
		}
	}

private:

	void unsafe_process() {
		process_gen_rows() ;
		process_sample_rows() ;
		construct_plots() ;
	}

	void process_gen_rows() {
		m_context.logger() << "Processing SNPs...\n" ;
		Timer timer ;

		InternalStorageGenRow row ;
		std::size_t number_of_snps_processed = 0 ;
		for(
			row.set_number_of_samples( m_context.snp_data_source().number_of_samples() ) ;
			row.read_from_source( m_context.snp_data_source() ) ;
			row.set_number_of_samples( m_context.snp_data_source().number_of_samples() )
		) {
			preprocess_gen_row( row ) ;
			process_gen_row( row, ++number_of_snps_processed ) ;
			accumulate_per_column_amounts( row, m_per_column_amounts ) ;
			m_context.print_progress_if_necessary() ;
		}
		m_context.logger() << "\n" ;
		assert( m_context.snp_data_source().number_of_snps_read() == m_context.snp_data_source().total_number_of_snps() ) ;
		m_context.logger() << "Processed " << m_context.snp_data_source().total_number_of_snps() << " SNPs in "
			<< std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
		if( m_context.snp_filter().number_of_subconditions() > 0 ) {
			m_context.logger() << "(" << m_context.fltrd_in_snp_data_sink().number_of_snps_written() << " of " << m_context.snp_data_source().total_number_of_snps() << " SNPs passed the filter.)\n" ;
		}
	}
	
	void preprocess_gen_row( InternalStorageGenRow& row ) const {
		check_gen_row( row ) ;
		row.filter_out_samples_with_indices( m_context.indices_of_filtered_out_samples() ) ;
	}
	
	void check_gen_row( GenRow& row ) const {
		check_gen_row_has_correct_number_of_samples( row ) ;
	}

	void check_gen_row_has_correct_number_of_samples( GenRow& row ) const {
		if( row.number_of_samples() != m_context.sample_rows().size() ) {
			throw GenAndSampleFileMismatchException( row, row.number_of_samples(), m_context.sample_rows().size() ) ;
		}
	}

	void process_gen_row( GenRow const& row, std::size_t row_number ) {
		m_context.snp_statistics().process( row ) ;
		if( m_context.snp_filter().check_if_satisfied( m_context.snp_statistics() )) {
			row.write_to_sink( m_context.fltrd_in_snp_data_sink() ) ;
			output_gen_row_stats( m_context.snp_statistics(), row_number ) ;
		}
		else {
			row.write_to_sink( m_context.fltrd_out_snp_data_sink() ) ;
			do_snp_filter_diagnostics( m_context.snp_statistics(), row_number ) ;
		}
	}

	void output_gen_row_stats( GenotypeAssayStatistics const& row_statistics, std::size_t row_number ) {
		if( m_context.snp_statistics().size() > 0 ) {
			m_context.snp_stats_sink() << row_number ;
			for( std::size_t i = 0 ; i < row_statistics.size(); ++i ) {
				m_context.snp_stats_sink() << row_statistics.get_value< std::string >( row_statistics.get_statistic_name( i )) ;
			}
			m_context.snp_stats_sink() << statfile::end_row() ;
		}
	}

	void do_snp_filter_diagnostics( GenRowStatistics const& row_statistics, std::size_t const row_index ) {
		std::ostream& log = m_context.logger()[ "log" ] ;
		log
			<< "Filtered out snp #" << row_index << " (" << row_statistics.row().SNPID() << " " << row_statistics.row().RSID() << " " << row_statistics.row().SNP_position() << ")"
			<< " because it does not satisfy " ;
		for( std::size_t i = 0, failed_condition_count = 0; i < m_context.snp_filter().number_of_subconditions() ; ++i ) {
			if( !m_context.snp_filter().subcondition(i).check_if_satisfied( row_statistics )) {
				++m_context.snp_filter_failure_counts()[ i ] ;
				if( failed_condition_count > 0 ) {
					log << " or " ;
				}
				log << "\"" << m_context.snp_filter().subcondition(i) << "\"" ;
				++failed_condition_count ;
			}
		}
		log << ".\n" ;
	}

	void accumulate_per_column_amounts( GenRow& row, std::vector< GenotypeProportions >& per_column_amounts ) {
		// Keep totals for per-column stats.
		if( per_column_amounts.empty() ) {
			per_column_amounts.reserve( row.number_of_samples() ) ;
			std::copy( row.begin_genotype_proportions(), row.end_genotype_proportions(), std::back_inserter( per_column_amounts )) ;
		}
		else {
			assert( per_column_amounts.size() == row.number_of_samples() ) ;
			std::transform( per_column_amounts.begin(), per_column_amounts.end(),
			 				row.begin_genotype_proportions(),
			 				per_column_amounts.begin(),
							std::plus< GenotypeProportions >() ) ;
		}
	}

	void process_sample_rows() {
		m_context.logger() << "Processing samples...\n" ;
		Timer timer ;

		apply_sample_filter() ;
		assert( m_context.sample_rows().size() == m_per_column_amounts.size() ) ;

		for( std::size_t i = 0 ; i < m_per_column_amounts.size(); ++i ) {
			SampleRow& sample_row = m_context.sample_rows()[i] ;
			m_context.sample_statistics().process( sample_row, m_per_column_amounts[i], m_context.snp_data_source().total_number_of_snps() ) ;
			output_sample_stats( i + 1, m_context.sample_statistics() ) ;
			
			m_context.sample_statistics().add_to_sample_row( sample_row, "missing" ) ;
			m_context.sample_statistics().add_to_sample_row( sample_row, "heterozygosity" ) ;
			m_context.fltrd_in_sample_sink() << sample_row ;
		}

		m_context.logger() << "Processed " << m_context.snp_data_source().number_of_samples() << " samples in " << std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
		if( m_context.sample_filter().number_of_subconditions() > 0 ) {
			m_context.logger() << "(" << m_context.sample_rows().size() << " of " << m_context.snp_data_source().number_of_samples() << " samples passed the filter.)\n" ;
		}
		m_context.logger() << "\n" ;
	}

	void output_sample_stats( std::size_t index, GenotypeAssayStatistics const& stats ) {
		m_context.sample_stats_sink() << index ;
		for( std::size_t i = 0; i < m_context.sample_statistics().size(); ++i ) {
			m_context.sample_stats_sink() << stats.get_value< std::string >( stats.get_statistic_name( i )) ;
		}
		m_context.sample_stats_sink() << statfile::end_row() ;
	}

	void apply_sample_filter() {
		for( std::size_t pre_filter_i = 0, post_filter_i = 0; post_filter_i < m_context.sample_rows().size(); ++pre_filter_i ) {
			if( sample_row_is_filtered_out( pre_filter_i ) ) {
				m_context.fltrd_out_sample_sink().write( m_context.sample_rows()[ post_filter_i] ) ;
				do_sample_filter_diagnostics( m_context.sample_rows()[ post_filter_i ], pre_filter_i ) ;
				m_context.sample_rows().erase( m_context.sample_rows().begin() + post_filter_i ) ;
			}
			else {
				++post_filter_i ;
			}
		}
	}

	bool sample_row_is_filtered_out( std::size_t const sample_row_index ) {
		return std::binary_search( m_context.indices_of_filtered_out_samples().begin(), m_context.indices_of_filtered_out_samples().end(), sample_row_index ) ;
	}

	void do_sample_filter_diagnostics( SampleRow const& sample_row, std::size_t const sample_row_index ) {
		std::ostream& log = m_context.logger()["log"] ;
		log
			<< "Filtered out sample row " << sample_row_index << " (" << sample_row.ID1() << " " << sample_row.ID2() << ")"
			<< " because it does not satisfy " ;
		for( std::size_t i = 0, failed_condition_count = 0; i < m_context.sample_filter().number_of_subconditions() ; ++i ) {
			if( !m_context.sample_filter().subcondition(i).check_if_satisfied( sample_row )) {
				++m_context.sample_filter_failure_counts()[ i ] ;
				if( failed_condition_count > 0 ) {
					log << " or " ;
				}
				log << "\"" << m_context.sample_filter().subcondition(i) << "\"" ;
				++failed_condition_count ;
			}
		}
		log << ".\n" ;
	}

	void construct_plots() {
		// not implemented.
	}

private:
	QCToolContext& m_context ;

	std::vector< GenotypeProportions > m_per_column_amounts ;
} ;

#endif
