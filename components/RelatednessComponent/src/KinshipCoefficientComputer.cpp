
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/function.hpp>
#include <boost/timer/timer.hpp>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include "unistd.h"
#include "config/config.hpp"
#if HAVE_CBLAS
	#include "cblas.h"
#endif

#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "worker/Task.hpp"
#include "worker/FunctionTask.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "components/RelatednessComponent/KinshipCoefficientComputer.hpp"
#include "components/RelatednessComponent/PCAComputer.hpp"
#include "components/RelatednessComponent/mean_centre_genotypes.hpp"

#define USING_BOOST_THREADPOOL 0
#if USING_BOOST_THREADPOOL
#include "boost/threadpool.hpp"
#else
#include "metro/concurrency/threadpool.hpp"
#endif

// #define DEBUG_KINSHIP_COEFFICIENT_COMPUTER 2

namespace impl {
	std::ostream& operator<<( std::ostream& out, SampleBounds const& b ) {
		out << "[" << b.begin_sample_i << "-" << b.end_sample_i << ", " << b.begin_sample_j << "-" << b.end_sample_j << "]" ;
		return out ;
	}
	
	// Return a tiling of a dxd matrix into (roughly) B blocks.
	// Blocks should be roughly evenly sized.
	// This implementation makes off-diagonal blocks half the size of diagonal blocks.
	std::vector< SampleBounds > get_matrix_lower_diagonal_tiling( std::size_t d, std::size_t B ) {
		// Because off-diagonal blocks involve about twice as much computation,
		// we make them half the size vertically.
		// If there are (N ( N+1)/2 ) square blocks this makes
		// (N (N+1) ) - N = N^2 blocks in total.
		// So we want to choose N so that N^2 is about twice the number of threads.
		// Def
		std::size_t const N = B ;
		std::vector< SampleBounds > result ;
		double K = double( d ) / N ;
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
		std::cerr << "get_matrix_lower_diagonal_tiling: " << d << "x" << d << " matrix, "
			<< N << " blocks across diagonal, "
			<< N*(N-1) << " blocks in lower triangle.\n" ;
//#endif
		// Set up some initial tasks.
		// Do j (the column) as the outer loop
		// so that tasks go in column-major direction.
		for( std::size_t j = 0; j < N; ++j ) {
			for( std::size_t i = 0; i < N; ++i ) {
				if( i == j ) {
					SampleBounds bounds ;
					bounds.begin_sample_i = std::ceil(i*K) ;
					bounds.end_sample_i = std::min( std::size_t( std::ceil( (i+1)*K ) ), d ) ;
					bounds.begin_sample_j = bounds.begin_sample_i ;
					bounds.end_sample_j = bounds.end_sample_i ;
					result.push_back( bounds ) ;
				}
				else {
					// Do two vertical stripes
					SampleBounds bounds ;
					bounds.begin_sample_i = std::ceil(i*K) ;
					bounds.end_sample_i = std::min( std::size_t( std::ceil( (i+1)*K ) ), d ) ;
					bounds.begin_sample_j = std::ceil(j*K) ;
					bounds.end_sample_j = std::ceil((j+0.5)*K) ;
					result.push_back( bounds ) ;

					bounds.begin_sample_j = bounds.end_sample_j ;
					bounds.end_sample_j = std::min( std::size_t( std::ceil( (j+1)*K ) ), d ) ;
					result.push_back( bounds ) ;
				}
			}
		}
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
		std::cerr << "get_matrix_lower_diagonal_tiling: tiling has " << result.size() << " blocks in total.\n" ;
//#endif
		return result ;
	}


	std::vector< SampleBounds > get_matrix_lower_diagonal_vertical_strip_tiling( std::size_t d, std::size_t const strip_width ) {
		std::vector< SampleBounds > result ;
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
		std::cerr << "get_matrix_lower_diagonal_vertical_strip_tiling: " << d << "x" << d << " matrix" ;
//#endif
		// Set up some initial tasks.
		// Do j (the column) as the outer loop
		// so that tasks go in column-major direction.
		for( std::size_t i = 0; i < d; i += strip_width ) {
			SampleBounds bounds ;
			bounds.begin_sample_i = i ;
			bounds.end_sample_i = d ;
			bounds.begin_sample_j = i ;
			bounds.end_sample_j = std::min( i+strip_width, d ) ;
			result.push_back( bounds ) ;
		}
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
		std::cerr << "get_matrix_lower_diagonal_vertical_strip_tiling: tiling has " << result.size() << " blocks in total.\n" ;
//#endif
		return result ;
	}
	
	std::vector< SampleBounds > get_matrix_lower_diagonal_horizontal_strip_tiling( std::size_t d, std::size_t const strip_width ) {
		std::vector< SampleBounds > result ;
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
		std::cerr << "get_matrix_lower_diagonal_horizontal_strip_tiling: " << d << "x" << d << " matrix" ;
//#endif
		// Set up some initial tasks.
		// Do j (the column) as the outer loop
		// so that tasks go in column-major direction.
		for( std::size_t i = 0; i < d; i += strip_width ) {
			SampleBounds bounds ;
			bounds.begin_sample_i = i ;
			bounds.end_sample_i = i + strip_width ;
			bounds.begin_sample_j = 0 ;
			bounds.end_sample_j = std::min( i+strip_width, d ) ;
			result.push_back( bounds ) ;
		}
//#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
		std::cerr << "get_matrix_lower_diagonal_horizontal_strip_tiling: tiling has " << result.size() << " blocks in total.\n" ;
//#endif
		return result ;
	}

	NTaskDispatcher::NTaskDispatcher(
		worker::Worker* worker
	):
		m_worker( worker )
	{
	}
	
	void NTaskDispatcher::set_number_of_tasks( std::size_t number_of_tasks ) {
		wait_until_complete() ;
		// TODO: do this better using the constructor
		// I don't have the docs here so don't know the right call.
		m_tasks.clear() ;
		for( std::size_t i = 0; i < number_of_tasks; ++i ) {
			m_tasks.push_back( 0 ) ;
		}
	}

	std::size_t NTaskDispatcher::number_of_tasks() const { return m_tasks.size() ; }
	
	void NTaskDispatcher::submit_task( std::size_t task_i, boost::function< void() > task ) {
		assert( task_i < m_tasks.size() ) ;
		if( !m_tasks.is_null( task_i ) ) {
			throw genfile::BadArgumentError(
				"impl::NTaskDispatcher::submit_task()",
				( boost::format( "task_i=%d" ) % task_i ).str(),
				"Task is still running.  Call wait_until_complete() first."
			) ;
		}
		worker::Task::UniquePtr new_task( new worker::FunctionTask( task ) ) ;
		m_tasks.replace( task_i, new_task ) ;
		m_worker->tell_to_perform_task( m_tasks[ task_i ] ) ;
	}
	
	void NTaskDispatcher::wait_until_complete() {
		for( std::size_t i = 0; i < m_tasks.size(); ++i ) {
			if( !m_tasks.is_null( i ) ) {
				m_tasks[ i ].wait_until_complete() ;
				m_tasks.replace( i, 0 ) ;
			}
		}
	}

	NormaliseGenotypesAndComputeXXtFast::UniquePtr NormaliseGenotypesAndComputeXXtFast::create(
		worker::Worker* worker,
		std::size_t const number_of_snps_per_computation
	) {
		return NormaliseGenotypesAndComputeXXtFast::UniquePtr( new NormaliseGenotypesAndComputeXXtFast( worker, number_of_snps_per_computation ) ) ;
	}

	NormaliseGenotypesAndComputeXXtFast::NormaliseGenotypesAndComputeXXtFast(
		worker::Worker* worker,
		std::size_t const number_of_snps_per_chunk
	):
		m_worker( worker ),
		m_dispatcher( new NTaskDispatcher( worker ) ),
		m_call_threshhold( 0.9 ),
		m_allele_frequency_threshhold( 0.001 ),
		m_number_of_snps_per_chunk( number_of_snps_per_chunk ),
		m_number_of_snps_per_computation( 4 * number_of_snps_per_chunk )
	{
		// We allocate lookup tables exactly once, here.
		// We follow, roughly, the scheme in plink 1.9 of having four lookup tables.
		// Here we encode pairs of genotypes in 4 bits (though we could do 3 as in plink).
		// Each computation then involves 16 variants, four from each lookup table.
		// This multi-table version reduces the number of times the result matrices
		// have to be traversed, and seems to make most difference in speed when
		// multithreading.
		m_lookup_tables.resize( 4, std::vector< double >( 1 << ( m_number_of_snps_per_chunk * 4 ), 0.0 ) ) ;
		m_nonmissingness_lookup_tables.resize( 4, std::vector< int >( 1 << ( m_number_of_snps_per_chunk * 4 ), 0.0 ) ) ;
	}

	NormaliseGenotypesAndComputeXXtFast::Matrix const& NormaliseGenotypesAndComputeXXtFast::result() const {
		return m_result ;
	}

	NormaliseGenotypesAndComputeXXtFast::IntegerMatrix const& NormaliseGenotypesAndComputeXXtFast::nonmissingness() const {
		return m_nonmissingness ;
	}

	std::size_t NormaliseGenotypesAndComputeXXtFast::number_of_snps_included() const {
		return m_snp_count ;
	}

	std::string NormaliseGenotypesAndComputeXXtFast::get_summary() const {
		return (
			boost::format(
				"NormaliseGenotypesAndComputeXXtFast():\n"
				" - %d SNPs per chunk\n"
				" - %d chunks per array acces\n"
				" - four lookup tables of size: %d\n"
				" - minimum allele frequency: %.3f"
			)
				% m_number_of_snps_per_chunk
				% (m_number_of_snps_per_computation / m_number_of_snps_per_chunk)
				% m_lookup_tables[0].size()
				% m_allele_frequency_threshhold
		).str()
		;
	}
	void NormaliseGenotypesAndComputeXXtFast::begin_processing_snps(
		std::size_t number_of_samples,
		genfile::SNPDataSource::Metadata const&
	) {
		m_result.resize( number_of_samples, number_of_samples ) ;
		m_result.setZero() ;
		m_nonmissingness.resize( number_of_samples, number_of_samples ) ;
		m_nonmissingness.setZero() ;
		std::size_t numberOfTasks = m_worker->get_number_of_worker_threads() ;
#if 0
		m_matrix_tiling = get_matrix_lower_diagonal_vertical_strip_tiling( number_of_samples, 1 ) ;
#else
		m_matrix_tiling = get_matrix_lower_diagonal_tiling(
			number_of_samples,
			numberOfTasks
		) ;
#endif
		m_dispatcher->set_number_of_tasks( m_matrix_tiling.size() ) ;
		m_combined_genotypes.resize( number_of_samples, 0 ) ;
		m_per_snp_genotypes.resize( number_of_samples, 0 ) ;
		m_snp_count = 0 ;
	}

	namespace impl {
		template< typename T >
		double add( T a, T b ) {
			return a + b ;
		}
	}

#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
	namespace {
		template< typename T >
		void print_lookup_table( std::ostream& out, std::vector< T > const& lookup_table, std::size_t const max = 64 ) {
			out << "print_lookup_table(): lookup table of size " << lookup_table.size() << ":" ;
			for( std::size_t i = 0; i < std::min( max, lookup_table.size() ); ++i ) {
				if( i % 8 == 0 ) {
					out << "\n" ;
				} else {
					out << " " ;
				}
				out << ( boost::format( " - 0x%02x: %.2f" ) % i % lookup_table[i] ) ;
			}
			out << "\n" ;
		}
	}
#endif

	void NormaliseGenotypesAndComputeXXtFast::processed_snp(
		genfile::VariantIdentifyingData const& id_data,
		genfile::VariantDataReader::SharedPtr data_reader
	) {
		// std::vector< std::size_t > m_genotypes
		// std::vector< int >
		m_per_snp_genotypes.assign( m_result.rows(), 0 ) ;

		// I find it simplest here to encode genotypes as
		// 0 (missing), 1 (AA homozygote), 2 (heterozygote), 3 (BB homozygote).
		data_reader->get(
			":genotypes:",
			genfile::vcf::get_threshholded_calls( m_per_snp_genotypes, m_call_threshhold, 0, 1, 2, 3 )
		) ;

		// compute allele frequency
		double allele2_count = 0.0 ;
		double nonmissing_count = 0.0 ;
		for( std::size_t i = 0; i < m_combined_genotypes.size(); ++i ) {
			bool const missing = ( m_per_snp_genotypes[i] == 0 ) ;
			if( !missing ) {
				allele2_count += m_per_snp_genotypes[i] - 1.0 ;
				nonmissing_count += 1 ;
			}
		}
		// Data estimate:
		// double const allele_frequency = allele2_count / ( nonmissing_count * 2.0 ) ;
		// Regularised estimate based on adding one of each allele,
		// this matches Price et al Nature Genetics 2006
		double const posterior_allele_frequency = (1 + allele2_count) / ( 2 + ( 2 * nonmissing_count)) ;
		double const allele_frequency = allele2_count / ( 2 * nonmissing_count ) ;
		if( nonmissing_count > 0 && std::min( allele_frequency, 1.0 - allele_frequency ) > m_allele_frequency_threshhold ) {
			// Ok we will process this SNP.
			// We batch SNPs together to form a total of m_number_of_snps_per_computation SNPs.
			// This means we visit the result matrix (which for large samples won't fit in cache)
			// as few times as possible.
			// Here we are not quite as efficient as plink: I use four bits instead of three to
			// encode the two genotypes.
			std::size_t const lookup_snp_index = m_snp_count % m_number_of_snps_per_computation ;
			if( lookup_snp_index == 0 ) {
				m_dispatcher->wait_until_complete() ;
				std::copy( m_per_snp_genotypes.begin(), m_per_snp_genotypes.end(), m_combined_genotypes.begin() ) ;
				// Clear the lookup tables
				for( std::size_t i = 0; i < m_lookup_tables.size(); ++i ) {
					std::fill( m_lookup_tables[i].begin(), m_lookup_tables[i].end(), 0.0 ) ;
					std::fill( m_nonmissingness_lookup_tables[i].begin(), m_nonmissingness_lookup_tables[i].end(), 0.0 ) ;
				}
			} else {
				// We allow 4 bits per genotype, although each genotype takes up only the lower two bits.
				// This allows us to encode pairs of genotypes in consecutive pairs of bits
				// for simpler lookup of the lookup table.
				for( std::size_t i = 0; i < m_combined_genotypes.size(); ++i ) {
					m_combined_genotypes[i] |= m_per_snp_genotypes[i] << (4 * lookup_snp_index) ;
				}
			}

			std::size_t const lookup_table_index = ( m_snp_count % m_number_of_snps_per_computation ) / m_number_of_snps_per_chunk ;
			double const mean = 2.0 * allele_frequency ;
			double const sd = std::sqrt( ( 2.0 * posterior_allele_frequency * ( 1.0 - posterior_allele_frequency )) ) ;
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
			std::cerr << "NormaliseGenotypesAndComputeXXtFast::processed_snp(): filling lookup table " << lookup_table_index << ".\n" ;
#endif
			
			add_snp_to_lookup_table(
				lookup_snp_index % m_number_of_snps_per_chunk,
				mean, sd,
				&m_lookup_tables[ lookup_table_index ],
				&m_nonmissingness_lookup_tables[ lookup_table_index ]
			) ;

			++m_snp_count ;
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
			boost::format fmt( "%.2f" ) ;
			std::cerr << "NormaliseGenotypesAndComputeXXtFast::processed_snp(): m_snp_count: " << m_snp_count << ".\n" ;
			std::cerr << "NormaliseGenotypesAndComputeXXtFast::processed_snp(): allele_frequency: " << fmt % allele_frequency << ".\n" ;
			std::cerr << "NormaliseGenotypesAndComputeXXtFast::processed_snp(): mean: " << fmt % mean << ", sd = " << fmt % sd << ".\n" ;
			print_lookup_table( std::cerr, m_lookup_tables[0], 32 ) ;
			print_lookup_table( std::cerr, m_lookup_tables[1], 32 ) ;
			print_lookup_table( std::cerr, m_lookup_tables[2], 32 ) ;
			print_lookup_table( std::cerr, m_lookup_tables[3], 32 ) ;
			//print_lookup_table( std::cerr, m_nonmissingness_lookup_table, 256 ) ;
#endif
			if( m_snp_count % m_number_of_snps_per_computation == 0 ) {
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER
			std::cerr << "NormaliseGenotypesAndComputeXXtFast::processed_snp(): submitting tasks...\n" ;
#endif
				submit_tasks( m_combined_genotypes ) ;
			}
		}
	}

	void NormaliseGenotypesAndComputeXXtFast::add_snp_to_lookup_table(
		std::size_t const snp_lookup_index,
		double const mean,
		double const sd,
		std::vector< double >* lookup_table,
		std::vector< int >* nonmissingness_lookup_table
	) {
		// We lay out the lookup table as follows.
		// The two genotypes at the first SNP occupy the first 4 bits,
		// the two genotypes at the second SNP occupy the next 4 bits,
		// etc.
		// Thus after k SNPs the size of the useful part of the table is 16^k,
		// which is computed as (1 << 4k).
		// Note that because genotype value 0 corresponds to missing data,
		// The block as computed after SNP k-1 is equal to the block
		// for missing values for SNP k.
		if( snp_lookup_index == 0 ) {
			// bottom 2 bits correspond to top genotype = 0.
			for( std::size_t g = 0; g < 4; ++g ) {
				(*nonmissingness_lookup_table)[g] = 0 ;
				(*lookup_table)[ g ] = 0 ;
			}
			for( std::size_t g = 1; g < 4; ++g ) {
				(*nonmissingness_lookup_table)[ (g << 2) + 0 ] = 0 ;
				(*nonmissingness_lookup_table)[ (g << 2) + 1 ] = 1 ;
				(*nonmissingness_lookup_table)[ (g << 2) + 2 ] = 1 ;
				(*nonmissingness_lookup_table)[ (g << 2) + 3 ] = 1 ;
				(*lookup_table)[ (g << 2) + 0 ] = 0.0 ;
				(*lookup_table)[ (g << 2) + 1 ] = ( ( 0.0 - mean ) / sd ) * ( ( (g-1) - mean ) / sd ) ;
				(*lookup_table)[ (g << 2) + 2 ] = ( ( 1.0 - mean ) / sd ) * ( ( (g-1) - mean ) / sd ) ;
				(*lookup_table)[ (g << 2) + 3 ] = ( ( 2.0 - mean ) / sd ) * ( ( (g-1) - mean ) / sd ) ;
			}
		} else {
			// We have already computed the initial block of size 16^k.
			// Now we compute 15 additional blocks for 16 in total.
			// In the computation we break these up as 4 blocks of size 
			// This means this SNP has no effect and we need only copy.
			std::size_t const inner_block_size = 1 << ( 4 * snp_lookup_index ) ;
			std::size_t const outer_block_size = 4 * inner_block_size ;
			for( std::size_t outer_g = 0; outer_g < 4; ++outer_g ) {
				for( std::size_t inner_g = (outer_g == 0) ? 1 : 0; inner_g < 4; ++inner_g ) {
					double const value = ( inner_g == 0 || outer_g == 0 )
						? 0.0
						: (((inner_g-1.0) - mean)/sd) * (((outer_g-1.0) - mean)/sd) ;
					int const nonmissingness = ( inner_g == 0 || outer_g == 0 )
						? 0
						: 1 ;
					std::transform(
						lookup_table->begin(),
						lookup_table->begin() + inner_block_size,
						lookup_table->begin() + outer_block_size * outer_g + inner_block_size * inner_g,
						boost::bind( &impl::add< double >, _1, value )
					) ;
					std::transform(
						nonmissingness_lookup_table->begin(),
						nonmissingness_lookup_table->begin() + inner_block_size,
						nonmissingness_lookup_table->begin() + outer_block_size * outer_g + inner_block_size * inner_g,
						boost::bind( &impl::add< int >, _1, nonmissingness )
					) ;
				}
			}
		}
	}

	void NormaliseGenotypesAndComputeXXtFast::submit_tasks(
		std::vector< uint64_t > const& genotypes
	) {
		for( std::size_t tile_i = 0; tile_i < m_matrix_tiling.size(); ++tile_i ) {
	#if 1
			m_dispatcher->submit_task(
				tile_i,
				boost::bind(
					&NormaliseGenotypesAndComputeXXtFast::compute_block,
					this,
					boost::cref( m_matrix_tiling[ tile_i ] ),
					boost::cref( genotypes )
				)
			) ;
	#else
			compute_block( m_matrix_tiling[tile_i], genotypes, lookup_table, nonmissingness_lookup_table ) ;
	#endif
		}
	}

	void NormaliseGenotypesAndComputeXXtFast::end_processing_snps() {
		if( m_snp_count % m_number_of_snps_per_computation > 0 ) {
			submit_tasks( m_combined_genotypes ) ;
		}
		m_dispatcher->wait_until_complete() ;
		std::cerr << "snp count is: " << m_snp_count << ".\n" ;
		std::cerr << "non-missingness array is:\n" << m_nonmissingness.block(0,0,10,10) << ".\n" ;
		m_result.array() /= m_nonmissingness.array().cast< double >() ;
	}

	void NormaliseGenotypesAndComputeXXtFast::compute_block(
		SampleBounds const& sample_bounds,
		std::vector< uint64_t > const& genotypes
	) {
		// only compute the lower diagonal.
		// We also want to make the function go down the rows of the result
		// in the inner loop, since it's stored column-major, so put row index in the inner loop.
		for( std::size_t j = sample_bounds.begin_sample_j; j < sample_bounds.end_sample_j; ++j ) {
			for( std::size_t i = std::max( j, std::size_t( sample_bounds.begin_sample_i ) ); i < sample_bounds.end_sample_i; ++i ) {
				uint64_t const combined_genotype = ( genotypes[i] << 2 ) + genotypes[j] ;
				std::size_t const lookup_index0 = combined_genotype & 0xFFFF ;
				std::size_t const lookup_index1 = ( combined_genotype >> 16 ) & 0xFFFF ;
				std::size_t const lookup_index2 = ( combined_genotype >> 32 ) & 0xFFFF ;
				std::size_t const lookup_index3 = ( combined_genotype >> 48 ) & 0xFFFF ;
				m_result(i,j) += m_lookup_tables[0][ lookup_index0 ] ;
				m_result(i,j) += m_lookup_tables[1][ lookup_index1 ] ;
				m_result(i,j) += m_lookup_tables[2][ lookup_index2 ] ;
				m_result(i,j) += m_lookup_tables[3][ lookup_index3 ] ;
				m_nonmissingness(i,j) += m_nonmissingness_lookup_tables[0][ lookup_index0 ] ;
				m_nonmissingness(i,j) += m_nonmissingness_lookup_tables[1][ lookup_index1 ] ;
				m_nonmissingness(i,j) += m_nonmissingness_lookup_tables[2][ lookup_index2 ] ;
				m_nonmissingness(i,j) += m_nonmissingness_lookup_tables[3][ lookup_index3 ] ;
			}
		}
	}
}

KinshipCoefficientComputer::KinshipCoefficientComputer(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	appcontext::UIContext& ui_context,
	Computation::UniquePtr computation
):
	m_options( options ),
	m_ui_context( ui_context ),
	m_samples( samples ),
	m_computation( computation )
{
	assert( m_computation.get() ) ;
}

void KinshipCoefficientComputer::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& metadata ) {
	m_computation->begin_processing_snps( number_of_samples, metadata ) ;
}

void KinshipCoefficientComputer::processed_snp( genfile::VariantIdentifyingData const& id_data, genfile::VariantDataReader::SharedPtr data_reader ) {
	m_computation->processed_snp( id_data, data_reader ) ;
}

void KinshipCoefficientComputer::end_processing_snps() {
	m_computation->end_processing_snps() ;

	std::string description = "Number of SNPs: "
		+ genfile::string_utils::to_string( m_computation->number_of_snps_included() )
		+ "\nNumber of samples: "
		+ genfile::string_utils::to_string( m_samples.get_number_of_individuals() )
	;
#if 1
	send_results(
		m_computation->number_of_snps_included(),
		m_computation->nonmissingness().cast< double >(),
		m_computation->result().cast< double >(),
		"KinshipCoefficientComputer",
		description
	) ;
	#endif
}
