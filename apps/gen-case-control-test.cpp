/*
 * This program, gen-case-control-test
 */

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <numeric>
#include <boost/math/special_functions/beta.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/bind.hpp>

#include <algorithm>

#include "ExternalStorageGenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "OptionProcessor.hpp"
#include "Whitespace.hpp"
#include "FileUtil.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SNPDataSource.hpp"
#include "SNPDataSourceChain.hpp"
#include "SNPDataSink.hpp"
#include "SNPDataSinkChain.hpp"
#include "string_utils.hpp"
#include "wildcard.hpp"
#include "parse_utils.hpp"
#include "Timer.hpp"
#include "progress_bar.hpp"

#include "CaseControlProbabilityPermuter.hpp"
#include "CaseControlPermutationsFileReader.hpp"
#include "CaseControlStatistic.hpp"
#include "BayesFactor.hpp"

namespace globals {
	std::string const program_name = "gen-case-control-test" ;
}

struct GenCaseControlProcessorException: public GToolException
{
	char const* what() const throw() { return "GenCaseControlProcessorException" ; }
} ;

struct GenCaseControlUsageError: public GenCaseControlProcessorException
{
	char const* what() const throw() { return "GenCaseControlUsageError" ; }
} ;

// thrown to indicate that the numbers of files specified for -cases or -controls
// is not two
struct GenCaseControlFileCountError: public GenCaseControlProcessorException
{
	char const* what() const throw() { return "GenCaseControlFileCountError" ; }
} ;

// thrown to indicate that an output file with wildcard appeared, but the corresponding input
// file had no wildcard.
struct GenCaseControlFileWildcardMismatchError: public GenCaseControlProcessorException
{
	char const* what() const throw() { return "GenCaseControlFileWildcardMismatchError" ; }
} ;	


struct GenCaseControlProcessor
{
public:
	static void declare_options( OptionProcessor & options ) {
		
		// File options		
	    options[ "-cases" ]
	        .set_description( "Path of gen file to input" )
			.set_is_required()
			.set_takes_values()
			.set_number_of_values_per_use( 1 ) ;

	    options[ "-controls" ]
	        .set_description( "Path of gen file to input" )
			.set_is_required()
			.set_takes_values()
			.set_number_of_values_per_use( 1 ) ;

		options [ "-permutations" ]
			.set_description( "Path of sample permutation file to input.")
			.set_is_required()
			.set_takes_single_value() ;

		options ["-use-differential-permuter"]
			.set_description( "Use the differential genotype probability permuter (which should be faster)")
			.set_default_value( false ) ;

		options [ "-force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
	}

	GenCaseControlProcessor( OptionProcessor const& options )
	:   m_options( options ),
		m_cout( std::cout.rdbuf() ),
	  	m_rng(0)
	{
		write_start_banner( m_cout ) ;
		setup() ;
	}

	~GenCaseControlProcessor() {
		write_end_banner( m_cout ) ;
	}
	
private:
	void setup() {
		try {
			get_required_filenames() ;
			open_gen_files() ;
			open_permutations_file() ;
			construct_rng() ;
		}
		catch( genfile::FileContainsSNPsOfDifferentSizes const& ) {
			m_cout << "The GEN files specified did not all have the same sample size.\n" ;
			throw ;
		}
		catch ( genfile::FileHasTwoConsecutiveNewlinesError const& e ) {
			m_cout << "\n!!ERROR: a GEN file was specified having two consecutive newlines.\n"
				<< "!! NOTE: popular editors, such as vim and nano, automatically add an extra newline to the file (which you can't see).\n"
				<< "!!     : Please check that each SNP in the file is terminated by a single newline.\n" ;
			throw ;
		}
	}

	void get_required_filenames() {
		std::vector< std::string > case_filenames, control_filenames ;

		get_gen_filenames( &m_control_gen_filenames, m_options.get_values< std::string >( "-controls" )) ;
		get_gen_filenames( &m_case_gen_filenames, m_options.get_values< std::string >( "-cases" )) ;
		m_permutations_filename = m_options.get_value< std::string >( "-permutations" ) ;
		assert( m_control_gen_filenames.size() > 0 ) ;
		assert( m_case_gen_filenames.size() > 0 ) ;
	}

	void get_gen_filenames( std::vector< std::string >* filename_list_ptr, std::vector< std::string > const& filenames ) {
		assert( filenames.size() > 0 ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			expand_and_add_filename( filename_list_ptr, filenames[i] ) ;
		}
	}

	void expand_and_add_filename( std::vector< std::string >* filename_list_ptr, std::string const& filename ) {
		bool input_file_has_wildcard = ( filename.find( '#' ) != std::string::npos ) ;
		if( input_file_has_wildcard ) {
			std::vector< wildcard::FilenameMatch >
				matches = wildcard::find_files_matching_path_with_integer_wildcard( filename, '#', 1, 100 ) ;
			for( std::size_t i = 0; i < matches.size(); ++i ) {
				add_filename( filename_list_ptr, matches[i].filename() ) ;
			}
		}
		else {
			add_filename( filename_list_ptr, filename ) ;
		}
	}

	bool check_if_string_is_a_number_from_1_to_100( std::string const& a_string ) {
		return parse_integer_in_half_open_range( a_string, 1, 101 ) != 101 ;
	}

	void add_filename( std::vector< std::string >* filename_list_ptr, std::string const& filename ) {
		filename_list_ptr->push_back( filename ) ;
	}

	void open_gen_files() {
		Timer timer ;
		open_gen_files( m_control_gen_filenames, m_control_gen_input_chain ) ;
		open_gen_files( m_case_gen_filenames, m_case_gen_input_chain ) ;

		if( timer.elapsed() > 1.0 ) {
			m_cout << "Opened " << m_control_gen_filenames.size() + m_case_gen_filenames.size() << " GEN files in " << timer.elapsed() << "s.\n" ;\
		}
	}
	
	void open_gen_files( std::vector< std::string > const& filenames, std::auto_ptr< genfile::SNPDataSourceChain >& chain ) {
		chain.reset( new genfile::SNPDataSourceChain() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			add_gen_file_to_chain( *chain, filenames[i] ) ;
		}
	}

	void add_gen_file_to_chain( genfile::SNPDataSourceChain& chain, std::string const& filename ) {
		Timer timer ;
		m_cout << "(Opening gen file \"" << filename << "\"...)" << std::flush ;
		std::auto_ptr< genfile::SNPDataSource > snp_data_source( genfile::SNPDataSource::create( filename )) ;
		chain.add_source( snp_data_source ) ;
		m_cout << " (" << timer.elapsed() << "s)\n" ;
	}

	void open_permutations_file() {
		Timer timer ;
		m_cout << "(Opening permutations file \"" << m_permutations_filename << "\"...)" << std::flush ;
		if( m_options.check_if_option_was_supplied( "-use-differential-permuter" )) {
			m_case_control_probability_permuter.reset( new DifferentialCaseControlProbabilityPermuter( m_permutations_filename )) ;
		}
		else {
			m_case_control_probability_permuter.reset( new SimpleCaseControlProbabilityPermuter( m_permutations_filename )) ;
		}
		m_cout << " (" << timer.elapsed() << "s)\n" ;
	}

	void construct_rng() {
		m_rng.reset( new boost::mt19937 ) ;
		m_cout << "mt19937: " << m_rng->min() << " " << m_rng->max() << ".\n" ;
	}
	
public:
	
	void write_start_banner( std::ostream& oStream ) const {
		oStream << "\nWelcome to " << globals::program_name << "\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}
	
	void write_end_banner( std::ostream& oStream ) const {
		oStream << "\n"
			<< "Thank you for using " << globals::program_name << ".\n" ;
	}

	void write_preamble( std::ostream& oStream ) const {
		oStream << std::string( 72, '=' ) << "\n\n" ;
		try {
			oStream << "Control GEN files:\n" ;
			print_gen_files( oStream, m_control_gen_filenames, m_control_gen_input_chain ) ;
			oStream << "Case GEN files:\n" ;
			print_gen_files( oStream, m_case_gen_filenames, m_case_gen_input_chain ) ;
			oStream << "Number of permutations: " << m_case_control_probability_permuter->number_of_permutations() << ".\n" ;
			oStream << std::string( 72, '=' ) << "\n\n" ;
		}
		catch (...) {
			oStream << std::string( 72, '=' ) << "\n\n" ;
			throw ;
		}
	}

	void print_gen_files( std::ostream& oStream, std::vector<std::string> const& filenames, std::auto_ptr< genfile::SNPDataSourceChain > const& chain ) const {
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			oStream
				<< "  (" << std::setw(6) << chain->number_of_snps_in_source( i ) << " snps)  "
				<< "\"" << std::setw(20) << filenames[i] << "\"\n" ;
		}
		oStream << "  (total " << chain->total_number_of_snps() << " snps).\n\n" ;
	}

	void do_checks() {
	}
	
	void process() {
		unsafe_process() ;
	}

private:

	void unsafe_process() {
		m_timer.restart() ;

		std::size_t number_of_control_samples = m_control_gen_input_chain->number_of_samples(),
			number_of_case_samples = m_case_gen_input_chain->number_of_samples() ;

		std::vector< GenotypeProbabilities > probabilities( number_of_control_samples + number_of_case_samples ) ;
		ExternalStorageGenRow
			control_row( &probabilities[0], number_of_control_samples ),
			case_row( &probabilities[ number_of_control_samples ], number_of_case_samples ) ;
		
		m_cout << "Processing SNPs...\n" ;

		SimpleBayesFactor< CachingBetaCalculator > bayes_factor ;
		m_last_time = -5.0 ;
		std::size_t number_of_control_snps_matched = 0 ;
		while( read_snp( *m_case_gen_input_chain, case_row )) {
			std::size_t number_of_matching_snps = 0 ;
			while( read_next_snp_with_specified_position( *m_control_gen_input_chain, control_row, case_row.SNP_position())) {
				++number_of_matching_snps ;
				++number_of_control_snps_matched ;
				process_snp( case_row, control_row, bayes_factor, probabilities ) ;

				print_progress_if_needed() ;
			}
			if( number_of_matching_snps != 1 ) {
				m_cout << "\nA strange number (" << number_of_matching_snps << ") of control snps matched the case snp "
					<< case_row.SNPID() << " "
					<< case_row.RSID() << " "
					<< case_row.SNP_position()
					<< ".\n" ;
			}
		}

		std::cerr
			<< "\nProcessed case / control data ("
			<< m_case_gen_input_chain->number_of_snps_read() << " case SNPs, "
			<< number_of_control_snps_matched << " matched control SNPs) in "
			<< std::fixed << std::setprecision(1) << m_timer.elapsed() << " seconds.\n" ;
	
		close_all_files() ;
	}
	
	void print_progress_if_needed() {
		double time_now = m_timer.elapsed() ;
		if( (time_now - m_last_time >= 1.0) || (m_case_gen_input_chain->number_of_snps_read() == m_case_gen_input_chain->total_number_of_snps()) ) {
			print_progress( time_now ) ;
			m_last_time = time_now ;
		}
	}
	
	void process_snp( GenRow const& case_row, GenRow const& control_row, CaseControlStatistic const& bayes_factor, std::vector< GenotypeProportions >& probabilities ) {
		// I think we'd expect all the SNP identifying data to be the same for the two snps.
		// If there are mismatches, print out some information about them.
		compare_snps( case_row, control_row ) ;
		calculate_bayes_factors( case_row, control_row, bayes_factor, probabilities ) ;
	}
	
	void compare_snps( GenRow const& case_row, GenRow const& control_row ) {
		if( case_row.SNPID() != control_row.SNPID() ) {
			m_cout << "\nMatching rows have differing SNPIDs " << case_row.SNPID() << " and " << control_row.SNPID() << ".\n" ;
		}
		if( case_row.RSID() != control_row.RSID() ) {
			m_cout << "\nMatching rows have differing RSIDs " << case_row.RSID() << " and " << control_row.RSID() << ".\n" ;
		}
		if( case_row.SNP_position() != control_row.SNP_position() ) {
			m_cout << "\nMatching rows have different positions " << case_row.SNP_position() << " and " << control_row.SNP_position() << ".\n" ;
		}
		if( case_row.first_allele() != control_row.first_allele() || case_row.second_allele() != control_row.second_allele() ) {
			m_cout 	<< "\nMatching rows have differing alleles "
					<< case_row.first_allele() << " " << case_row.second_allele()
					<< " and " << control_row.first_allele() << " " << control_row.second_allele() << ".\n" ;
		}
	}

	void calculate_bayes_factors( GenRow const& control_row, GenRow const& case_row, CaseControlStatistic const& bayes_factor, std::vector< GenotypeProbabilities >& probabilities )
	{
		std::map< double, GenotypeAmounts > status_to_genotype_amount_map ;

		double BF ;
		double average_BF = 0.0 ;
		Timer bayes_factor_timer ;

		GenotypeAssayStatistics stats1, stats2 ;

		std::size_t count = 0 ;

		m_case_control_probability_permuter->reset( probabilities ) ;
		do {
			print_progress_if_needed() ;

#if 0
			stats1.process( control_row.begin_genotype_proportions(), control_row.end_genotype_proportions() ) ;
			stats2.process( case_row.begin_genotype_proportions(), case_row.end_genotype_proportions() ) ;
			status_to_genotype_amount_map[0.0] = stats1.get_genotype_amounts() ; //m_case_control_probability_permuter->get_control_genotype_amounts() ;
			status_to_genotype_amount_map[1.0] = stats2.get_genotype_amounts() ; //m_case_control_probability_permuter->get_case_genotype_amounts() ;
#else
			status_to_genotype_amount_map[0.0] = m_case_control_probability_permuter->get_control_genotype_amounts() ;
			status_to_genotype_amount_map[1.0] = m_case_control_probability_permuter->get_case_genotype_amounts() ;
#endif
			double this_BF = bayes_factor.get_value< double >( status_to_genotype_amount_map ) ;
			if( count == 0 ) { //m_case_control_probability_permuter->current_permutation() == 0 ) {
				BF = this_BF ;
			}
			average_BF += this_BF ;
		}
//		while( ++count < m_case_control_probability_permuter->number_of_permutations() ) ; //m_case_control_probability_permuter->move_to_next_permutation() ) ;
		while( m_case_control_probability_permuter->move_to_next_permutation() ) ;

		m_cout << "\nComputed " << m_case_control_probability_permuter->number_of_permutations() << " bayes factors in "
			<< std::fixed << std::setprecision(1) << bayes_factor_timer.elapsed()
			<< "s (" << (m_case_control_probability_permuter->number_of_permutations() / bayes_factor_timer.elapsed()) << " BFs per sec).\n" ;
		average_BF /= (m_case_control_probability_permuter->number_of_permutations()) ;
		m_cout << "BF = " << BF << ", average BF = " << average_BF << ".\n" ;
	}

	void permute_probabilities( std::vector< GenotypeProbabilities >& probabilities ) {
		cyclically_permute_probabilities( probabilities ) ;
	}

	void randomly_permute_probabilities( std::vector< GenotypeProbabilities >& probabilities ) {
		// Do a random permutation, using Knuth's shuffle as described on the wikipedia page.
		std::size_t N = probabilities.size() ;
		GenotypeProbabilities* first_prob = &probabilities[0] ;
		for( std::size_t i = 1; i <= N; ++i ) {
			std::size_t k = get_random_number( N - i );
			if( k != ( N - i )) {
				std::swap( *(first_prob + k), *(first_prob + (probabilities.size() - i))) ;
			}
		}
	}

	// Get a random integer (hopefully) uniformly distributed in the range [0, L]
	std::size_t get_random_number( std::size_t L ) {
		// We take the output of the random number generator (which lies in the range [0, m_rng->max()])
		// Then multiply by L / m_rng->max().
		// Note that there are probably issues with the uniformity of the resulting number.
		return (static_cast< double >(L) / m_rng->max()) * static_cast< double >((*m_rng)()) ;
	}

	void cyclically_permute_probabilities( std::vector< GenotypeProbabilities >& probabilities ) {
		std::rotate( &probabilities[0], &probabilities[1], (&probabilities[0]) + probabilities.size() ) ;
	}

	void print_progress( double time_now ) {
		double case_progress = (static_cast< double >( m_case_gen_input_chain->number_of_snps_read() ) / m_case_gen_input_chain->total_number_of_snps()) ;
		double control_progress = (static_cast< double >( m_control_gen_input_chain->number_of_snps_read() ) / m_control_gen_input_chain->total_number_of_snps()) ;
		m_cout
			<< "\r"
			<< get_progress_bar( 30, case_progress )
			<< " (" << m_case_gen_input_chain->number_of_snps_read() << "/" << m_case_gen_input_chain->total_number_of_snps() << ")"
			<< " "
			<< get_progress_bar( 30, control_progress )
			<< " (" << m_control_gen_input_chain->number_of_snps_read() << "/" << m_control_gen_input_chain->total_number_of_snps() << ")"
			<< " (" << std::fixed << std::setprecision(1) << time_now << "s)"
			<< std::string( std::size_t(5), ' ' )
			<< std::flush ;
	}

	typedef boost::function< bool( std::string const&, std::string const&, uint32_t, char, char ) > SNPMatcher ;

	bool read_snp( genfile::SNPDataSource& snp_data_source, GenRow& row ) {
		return snp_data_source.read_snp(
			boost::bind< void >( &GenRow::set_number_of_samples, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNPID, &row, _1 ),
			boost::bind< void >( &GenRow::set_RSID, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNP_position, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele1, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele2, &row, _1 ),
			boost::bind< void >( &GenRow::set_genotype_probabilities, &row, _1, _2, _3, _4 )
		) ;
	}

	bool read_next_snp_with_specified_position( genfile::SNPDataSource& snp_data_source, GenRow& row, uint32_t specified_SNP_position ) {
		return snp_data_source.read_next_snp_with_specified_position(
			boost::bind< void >( &GenRow::set_number_of_samples, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNPID, &row, _1 ),
			boost::bind< void >( &GenRow::set_RSID, &row, _1 ),
			boost::bind< void >( &GenRow::set_SNP_position, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele1, &row, _1 ),
			boost::bind< void >( &GenRow::set_allele2, &row, _1 ),
			boost::bind< void >( &GenRow::set_genotype_probabilities, &row, _1, _2, _3, _4 ),
			specified_SNP_position 
		) ;
	}

	void close_all_files() {
		m_case_gen_input_chain.reset() ;
		m_control_gen_input_chain.reset() ;
	}
	
private:
	// Program options
	OptionProcessor const& m_options ;
	
	// filenames and files
	std::vector< std::string > m_case_gen_filenames ;
	std::vector< std::string > m_control_gen_filenames ;
	std::string m_permutations_filename ;
	std::auto_ptr< genfile::SNPDataSourceChain > m_case_gen_input_chain ;
	std::auto_ptr< genfile::SNPDataSourceChain > m_control_gen_input_chain ;
	std::ostream m_cout ;

	// Program operation
	std::vector< std::string > m_errors ;
	Timer m_timer ;
	double m_last_time ;
	std::auto_ptr< boost::mt19937 > m_rng ;

	// Genotype probability related
	std::vector< GenotypeProbabilities > m_probabilities ;
	std::auto_ptr< CaseControlProbabilityPermuter > m_case_control_probability_permuter ;
} ;


int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		GenCaseControlProcessor::declare_options( options ) ;
		options.process( argc, argv ) ;
    }
    catch( std::exception const& exception ) {
        std::cerr << "!! Error: " << exception.what() << ".\n";
        std::cerr << "Usage: " << globals::program_name << " [options]\n"
                << options
                << "\n" ;
        return -1 ;
    }

	// main section

	try	{
		GenCaseControlProcessor processor( options ) ;
		processor.do_checks() ;
		processor.write_preamble( std::cout ) ;
		processor.process() ;
	}
	catch( std::exception const& e )
	{
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}
