#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include "OptionProcessor.hpp"
#include "random_permutation.hpp"

namespace globals {
	std::string const program_name = "generate_random_case_control_status_file" ;
}

struct RandomCaseControlStatusGenerator
{
	static void declare_options( OptionProcessor& options ) {
		options[ "-number-of-cases" ]
			.set_takes_single_value()
			.set_is_required() ;

		options[ "-number-of-controls" ]
			.set_takes_single_value()
			.set_is_required() ;
			
		options[ "-number-of-permutations" ]
			.set_takes_single_value()
			.set_is_required()
			.set_default_value( 1000 ) ;
	}

	RandomCaseControlStatusGenerator( OptionProcessor const& options )
		: m_options( options )
	{
		write_start_banner( std::cerr ) ;
		setup() ;
	}

	void write_start_banner( std::ostream& oStream ) const {
		oStream << "\nWelcome to " << globals::program_name << "\n"
		 	<< "(C) 2009 University of Oxford\n\n";
	}
	
	~RandomCaseControlStatusGenerator()
	{
		write_end_banner( std::cerr ) ;
	}

	void write_end_banner( std::ostream& oStream ) const {
		oStream << "\n"
			<< "Thank you for using " << globals::program_name << ".\n" ;
	}

	void setup() {
		m_number_of_cases = m_options.get_value< std::size_t >( "-number-of-cases" ) ;
		m_number_of_controls = m_options.get_value< std::size_t >( "-number-of-controls" ) ;
		m_number_of_permutations = m_options.get_value< std::size_t >( "-number-of-permutations" ) ;
		assert( m_number_of_permutations > 1 ) ;
		m_case_control_status_vector.resize( m_number_of_cases + m_number_of_controls ) ;
	}

	void process() {
		std::cout << "Description: This file contains " << m_number_of_permutations << " random permutations of the 0-1 vector in the fourth line.\n" ;
		std::cout
			<< "number-of-zeroes: " << m_number_of_cases << "\n"
			<< "number-of-ones: " << m_number_of_controls << "\n"
			<< "number-of-permutations: " << m_number_of_permutations << "\n" ;
		construct_initial_state() ;
		print_state() ;
		for( std::size_t i = 0; i < m_number_of_permutations ; ++i ) {
			randomly_permute_range( m_case_control_status_vector.begin(), m_case_control_status_vector.end(), &generate_random_number_in_range_0_to_L ) ;
			print_state() ;
		}
	}

private:
	// Get a random integer (hopefully) uniformly distributed in the range [0, L]
	static std::size_t generate_random_number_in_range_0_to_L( std::size_t L ) {
		return (static_cast< double >(L) / m_rng.max()) * static_cast< double >(m_rng()) ;
	}

	void construct_initial_state() {
		std::size_t i = 0 ;
		for( ; i < m_number_of_controls; ++i ) {
			m_case_control_status_vector[i] = '0' ;
		}
		for( ; i < m_number_of_controls + m_number_of_cases; ++i ) {
			m_case_control_status_vector[i] = '1' ;
		}
	}

	void print_state() {
		for( std::size_t i = 0 ; i < (m_number_of_cases + m_number_of_controls) ; ++i ) {
			if( i > 0 ) {
				std::cout << ' ' ;
			}
			std::cout << m_case_control_status_vector[i] ;
		}
		std::cout << '\n' ;
	}

private:
	OptionProcessor const& m_options ;
	std::size_t m_number_of_cases, m_number_of_controls, m_number_of_permutations ;
	std::vector< char > m_case_control_status_vector ;
	static boost::mt19937 m_rng ;
} ;

boost::mt19937 RandomCaseControlStatusGenerator::m_rng ;

int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		RandomCaseControlStatusGenerator::declare_options( options ) ;
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
		RandomCaseControlStatusGenerator processor( options ) ;
		processor.process() ;
	}
	catch( std::exception const& e )
	{
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}

