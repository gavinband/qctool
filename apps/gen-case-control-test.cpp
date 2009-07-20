/*
 * This program, gen-case-control-test
 */

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <numeric>
#include <boost/math/special_functions/beta.hpp>

#include "ExternalStorageGenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "OptionProcessor.hpp"
#include "Whitespace.hpp"
#include "FileUtil.hpp"
#include "GenRowSource.hpp"
#include "GenRowSink.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "SNPDataSource.hpp"
#include "SNPDataSourceChain.hpp"
#include "SNPDataSink.hpp"
#include "SNPDataSinkChain.hpp"
#include "string_utils.hpp"
#include "wildcard.hpp"
#include "parse_utils.hpp"
#include "Timer.hpp"

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

std::vector< std::string > try_to_put_sample_file_first( std::string const& option_name, std::vector< std::string > const& filenames ) {
	std::vector< std::string > result( filenames ) ;
	if( filenames.size() != 2 ) {
		throw OptionValueInvalidException( option_name, filenames, "Two files, a sample and a gen file, must be specified for this option." ) ;
	}
	// swap the two options if it looks like the gen file is last.	
	if( !genfile::filename_indicates_gen_or_bgen_format( filenames[1] ) && genfile::filename_indicates_gen_or_bgen_format( filenames[0] )) {
		// Guess that sample file is listed second.
		std::swap( result[0], result[1] ) ;
	}
	return result ;
}


// Base class for individual statistics
struct CaseControlStatistic
{
	public:
		CaseControlStatistic() {} ;
		virtual ~CaseControlStatistic() {}

		typedef std::map< double, GenotypeAssayStatistics > case_status_statistic_map_t ;

		template< typename T >
		T get_value( case_status_statistic_map_t const& ) const ;

	protected:
		virtual double calculate_value( case_status_statistic_map_t const& ) const = 0;
} ;


template<>
double CaseControlStatistic::get_value< double >( case_status_statistic_map_t const& case_status_statistic_map ) const
{
	return calculate_value( case_status_statistic_map ) ;
}


//
// Calculate a bayes factor for case/control data sets,
// as on p.15 of  the 'Supplementary Methods' paper for
// "A new multipoint method for genome-wide association studies...", Marchini et al, 2007.
//
struct SimpleBayesFactor: public CaseControlStatistic
{
	SimpleBayesFactor(
		double phi_0 = 1.0,
		double eta_0 = 1.0,
		double phi_1 = 1.0,
		double eta_1 = 1.0,
		double phi_2 = 1.0,
		double eta_2 = 1.0
	)
		: m_phi(3), m_eta(3), m_precomputed_beta_of_phi_and_eta(3)
	{
		m_phi[0] = phi_0 ;
		m_eta[0] = eta_0 ;
		m_phi[1] = phi_1 ;
		m_eta[1] = eta_1 ;
		m_phi[2] = phi_2 ;
		m_eta[2] = eta_2 ;
		
		for( std::size_t i = 0; i < 3u; ++i ) {
			m_precomputed_beta_of_phi_and_eta[i] = beta( m_phi[i], m_eta[i] ) ;
		}
	}
	
protected:
	
	double calculate_value( case_status_statistic_map_t const& case_status_statistic_map ) const {
		// Assume we have map keys 0.0 and 1.0, and ignore all others.
		case_status_statistic_map_t::const_iterator
			control_i = case_status_statistic_map.find( 0.0 ),
			case_i = case_status_statistic_map.find( 1.0 ) ;

		GenotypeProbabilities const& control_genotype_amounts = control_i->second.get_genotype_amounts() ;
		GenotypeProbabilities const& case_genotype_amounts = case_i->second.get_genotype_amounts() ;
		return probability_of_data_given_M4( control_genotype_amounts, case_genotype_amounts )
			/ probability_of_data_given_M0( control_genotype_amounts, case_genotype_amounts ) ;
	}

	double beta( double a, double b ) const {
		return boost::math::beta( a, b ) ;
	}
	
	double probability_of_data_given_M4( GenotypeProbabilities const& control_genotype_amounts, GenotypeProbabilities const& case_genotype_amounts ) const {
		double result_for_AA_to_0_coding = 0.5;
		double result_for_BB_to_0_coding = 0.5 ;
		for( std::size_t i = 0; i < 3u; ++i ) {
			result_for_AA_to_0_coding *= beta( case_genotype_amounts[i] + m_phi[i], control_genotype_amounts[i] + m_eta[i] ) / m_precomputed_beta_of_phi_and_eta[i] ;
			result_for_BB_to_0_coding *= beta( case_genotype_amounts[2-i] + m_phi[i], control_genotype_amounts[2-i] + m_eta[i] ) / m_precomputed_beta_of_phi_and_eta[i] ;
		}
		return result_for_AA_to_0_coding + result_for_BB_to_0_coding ;
	}

	double probability_of_data_given_M0( GenotypeProbabilities const& control_genotype_amounts, GenotypeProbabilities const& case_genotype_amounts ) const {
		return beta( case_genotype_amounts.sum() + m_phi[0], control_genotype_amounts.sum() + m_eta[0] )
			/ m_precomputed_beta_of_phi_and_eta[0] ;
	}
	
private:
	
	std::vector< double > m_phi ;
	std::vector< double > m_eta ;
	std::vector< double > m_precomputed_beta_of_phi_and_eta ;
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
			.set_maximum_number_of_repeats( 1 )
			.set_number_of_values_per_use( 2 )
			.add_value_preprocessor( &try_to_put_sample_file_first ) ;

	    options[ "-controls" ]
	        .set_description( "Path of gen file to input" )
			.set_is_required()
			.set_takes_values()
			.set_maximum_number_of_repeats( 1 )
			.set_number_of_values_per_use( 2 )
			.add_value_preprocessor( &try_to_put_sample_file_first ) ;

		options[ "-number-of-permutations" ]
			.set_description( "The number of random permutations of case-control status to carry out at each SNP")
			.set_takes_single_value()
			.set_default_value( 1000 ) ;

		options [ "--force" ] 
			.set_description( "Ignore warnings and proceed with requested action." ) ;
	}

	GenCaseControlProcessor( OptionProcessor const& options )
		: m_cout( std::cout.rdbuf() ),
		  m_options( options )
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
			open_sample_files() ;
			// m_case_control_test.reset( CaseControlTest::create( m_control_sample_source, m_control_gen_input_chain, m_case_sample_source, m_case_gen_input_chain )) ;
			m_number_of_permutations = m_options.get_value<std::size_t >( "-number-of-permutations" ) ;
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

		control_filenames = m_options.get_values< std::string >( "-controls" ) ;
		assert( control_filenames.size() == 2 ) ;
		add_filename( &m_control_sample_filenames, control_filenames[0] ) ;
		expand_and_add_filename( &m_control_gen_filenames, control_filenames[1] ) ;

		case_filenames = m_options.get_values< std::string >( "-cases" ) ;
		assert( case_filenames.size() == 2 ) ;
		add_filename( &m_case_sample_filenames, case_filenames[0] ) ;
		expand_and_add_filename( &m_case_gen_filenames, case_filenames[1] ) ;
	}

	void expand_and_add_filename( std::vector< std::string >* filename_list_ptr, std::string const& filename ) {
		bool input_file_has_wildcard = ( filename.find( '#' ) != std::string::npos ) ;
		if( input_file_has_wildcard ) {
			std::pair< std::vector< std::string >, std::vector< std::string > >
				expanded_filename = find_files_matching_path_with_wildcard( filename, '#' ) ;

			// we only use matching filenames if the match is a number from 1 to 100
			// For such filenames, we add the filename to our list for cases.
			for( std::size_t j = 0; j < expanded_filename.first.size(); ++j ) {
				if( check_if_string_is_a_number_from_1_to_100( expanded_filename.second[j] )) {
					add_filename( filename_list_ptr, expanded_filename.first[j] ) ;
				}
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

	void open_sample_files() {
		assert( m_control_sample_filenames.size() == 1 ) ;
		Timer timer ;
		// m_control_sample_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_control_sample_filenames[0] ))) ;
		// m_case_sample_source.reset( new SampleInputFile< SimpleFileObjectSource< SampleRow > >( open_file_for_input( m_case_sample_filenames[0] ))) ;
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
		Timer timer ;

		std::size_t number_of_control_samples = m_control_gen_input_chain->number_of_samples(),
			number_of_case_samples = m_case_gen_input_chain->number_of_samples() ;

		std::vector< GenotypeProbabilities > probabilities( number_of_control_samples + number_of_case_samples ) ;
		ExternalStorageGenRow
			control_row( &probabilities[0], number_of_control_samples ),
			case_row( &probabilities[ number_of_control_samples ], number_of_case_samples ) ;
		
		m_cout << "Processing SNPs...\n" ;

		SimpleBayesFactor bayes_factor ;
		std::map< double, GenotypeAssayStatistics > case_control_statistic_map ;
		case_control_statistic_map[0.0] = GenotypeAssayStatistics() ;
		case_control_statistic_map[1.0] = GenotypeAssayStatistics() ;
		
		double last_time = -5.0 ;
		std::size_t number_of_control_snps_matched = 0 ;
		while( read_snp( *m_case_gen_input_chain, case_row )) {
			std::size_t number_of_matching_snps = 0 ;
			while( read_next_snp_with_specified_position( *m_control_gen_input_chain, control_row, case_row.SNP_position())) {
				++number_of_matching_snps ;
				++number_of_control_snps_matched ;
				// I think we'd expect all the SNP identifying data to be the same for the two snps.
				// If there are mismatches, print out some information about them.
				compare_snps( case_row, control_row ) ;

				case_control_statistic_map[0.0].process( control_row.begin_genotype_proportions(), control_row.end_genotype_proportions() ) ;
				case_control_statistic_map[1.0].process( case_row.begin_genotype_proportions(), case_row.end_genotype_proportions()	 ) ;

				double BF = bayes_factor.get_value< double >( case_control_statistic_map ) ;
				m_cout << "\nCalculated bayes factor: " << std::fixed << std::setprecision(6) << BF << ".\n" ;
				{
					Timer bayes_factor_timer ;
					for( std::size_t i = 0; i < m_number_of_permutations; ++i ) {
						permute_probabilities( probabilities ) ;
						case_control_statistic_map[0.0].process( control_row.begin_genotype_proportions(), control_row.end_genotype_proportions() ) ;
						case_control_statistic_map[1.0].process( case_row.begin_genotype_proportions(), case_row.end_genotype_proportions()	 ) ;
						
						bayes_factor.get_value< double >( case_control_statistic_map ) ;
					}	

					m_cout << "Computed " << m_number_of_permutations << " bayes factors in "
						<< std::fixed << std::setprecision(1) << bayes_factor_timer.elapsed()
						<< "s (" << (m_number_of_permutations / bayes_factor_timer.elapsed()) << " BFs per sec).\n" ;
				}

				// print a progress message every second.
				double time_now = timer.elapsed() ;
				if( (time_now - last_time >= 1.0) || (m_case_gen_input_chain->number_of_snps_read() == m_case_gen_input_chain->total_number_of_snps()) ) {
					print_progress( time_now ) ;
					last_time = time_now ;
				}
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
			<< std::fixed << std::setprecision(1) << timer.elapsed() << " seconds.\n" ;
	
		close_all_files() ;
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

	void permute_probabilities( std::vector< GenotypeProbabilities >& probabilities ) {
		// do nothing for now.
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

	std::string get_progress_bar( std::size_t width, double progress ) {
		progress = std::min( std::max( progress, 0.0 ), 1.0 ) ;
		std::size_t visible_progress = progress * width ;
		return
			"["
			+ std::string( std::size_t( visible_progress ), '*' )
			+ std::string( std::size_t( width - visible_progress ), ' ' )
			+ "]" ;
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
	
	std::ostream m_cout ;
	
	std::auto_ptr< genfile::SNPDataSourceChain > m_case_gen_input_chain ;
	std::auto_ptr< genfile::SNPDataSourceChain > m_control_gen_input_chain ;

	std::size_t m_number_of_permutations ;

	OptionProcessor const& m_options ;
	
	std::vector< std::string > m_case_gen_filenames ;
	std::vector< std::string > m_case_sample_filenames ;
	std::vector< std::string > m_control_gen_filenames ;
	std::vector< std::string > m_control_sample_filenames ;

	std::vector< std::string > m_errors ;
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
