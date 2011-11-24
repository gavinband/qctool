/*
 * This program, bechmark-statistics, reads some GEN rows from a file an
 *
 * Program arguments:
 *
 */

#include <numeric>
#include <string>
#include <sstream>
#include "../config.hpp"
#if HAVE_BOOST_TIMER
	#include <boost/timer.hpp>
	typedef boost::timer Timer ;
#else
	struct Timer
	{
		double elapsed() const { return 0.0 ;}
	} ;
#endif
#include "genfile/SNPDataSource.hpp"


// The following section defines the needed objects for use with the bgen.hpp implementation.
template< typename T >
struct Setter
{
	Setter( T& field ): m_field( field ) {} ;
	void operator()( T const& value ) { m_field = value ; }
private:
	T& m_field ;
} ;

template< typename T >
Setter< T > make_setter( T& field ) { return Setter<T>( field ) ; }


struct probabilities_t {
	double AA, AB, BB ;
	bool operator==( probabilities_t const& other ) const {
		return AA == other.AA && AB == other.AB && BB == other.BB ;
	}
} ;

struct ProbabilitySetter {
	ProbabilitySetter( std::vector< probabilities_t >& probabilities ): m_probabilities( probabilities ) {}
	void operator() ( std::size_t i, double aa, double ab, double bb ) {
		if( m_probabilities.size() < (i + 1) ) {
			m_probabilities.resize( i + 1) ;
		}
		m_probabilities[i].AA = aa ;  
		m_probabilities[i].AB = ab ;  
		m_probabilities[i].BB = bb ;  
	}

private:
	std::vector<probabilities_t>& m_probabilities ;
} ;

struct ProbabilityGetter {
	ProbabilityGetter( std::vector< probabilities_t > const& probabilities, int index ): m_probabilities( probabilities ), m_index( index ) {}
	double operator() ( std::size_t i ) const {
		switch( m_index ) {
			case 0: return m_probabilities[i].AA ; break ;
			case 1: return m_probabilities[i].AB ; break ;
			case 2: return m_probabilities[i].BB ; break ;
			default:
			assert(0); 
			break ;
		}
	}
private:
	std::vector<probabilities_t> const& m_probabilities ;
	int m_index ;
} ;

struct SnpData {
	
	SnpData() {} ;
	
	uint32_t number_of_samples ;
	std::string SNPID, RSID ;
	genfile::Chromosome chromosome ;
	uint32_t SNP_position ;
	std::string allele1, allele2 ;
	std::vector< probabilities_t > probabilities ;
	
	bool operator==( SnpData const& other ) const {
		return number_of_samples == other.number_of_samples
			&& SNPID == other.SNPID
			&& RSID == other.RSID
			&& chromosome == other.chromosome
			&& SNP_position == other.SNP_position
			&& allele1 == other.allele1
			&& allele2 == other.allele2
			&& probabilities == other.probabilities ;
	}
} ;

double process_gen_file( genfile::SNPDataSource& snp_data_source, std::size_t number_of_snps_to_read ) {
	SnpData snp_data ;
	std::cout << "Reading " << number_of_snps_to_read << " snps...\n" << std::flush ;
	std::size_t count = 0;
	Timer timer ;

	while( ( count < number_of_snps_to_read )
		&& snp_data_source.read_snp(
		make_setter( snp_data.number_of_samples ),
		make_setter( snp_data.SNPID ),
		make_setter( snp_data.RSID ),
		make_setter( snp_data.chromosome ),
		make_setter( snp_data.SNP_position ),
		make_setter( snp_data.allele1 ), 
		make_setter( snp_data.allele2 ),
		ProbabilitySetter( snp_data.probabilities )
	)) {
		++count ;
	}
	
	double elapsed = timer.elapsed() ;
	std::cout << "Read " << count << " snps in " << elapsed << "s.\n" ;
	return elapsed ;
}


int main( int argc, char** argv ) {
	if( argc < 2 ) {
		std::cerr << "You must supply an argument, the name of the gen file to load.\n" ;
		return -1 ;
	} else if( argc > 3 ) {
		std::cerr << "Too many arguments.\n" ;
		return -1 ;
	}

	// main section

	try	{
		std::string genFileName = argv[1] ;
		std::auto_ptr< genfile::SNPDataSource > source
			= genfile::SNPDataSource::create_chain( genfile::wildcard::find_files_by_chromosome( genFileName )) ;

		int number_of_snps_to_read = 10000 ;
		if( argc == 3 ) {
			std::istringstream stream ;
			stream.str( argv[2] ) ;
			stream >> number_of_snps_to_read ;
			if( !stream || number_of_snps_to_read <= 1 ) {
				std::cerr << "The second argument must be a positive integer.\n" ;
				return -1 ;
			}
		}

		std::size_t const number_of_tries = 3 ;
		std::vector< double > tries( number_of_tries ) ;
		for( std::size_t i = 0; i < number_of_tries; ++i ) {
			source->reset_to_start() ;
			tries[i] = process_gen_file( *source, number_of_snps_to_read ) ;
		}
		
		std::cerr << "\nAverage time taken was "
			<< std::accumulate( tries.begin(), tries.end(), 0.0 ) / number_of_tries
			<< ".\n" ;
	}
	catch( std::exception const& e )
	{
		std::cerr << "!! Error: " << e.what() << ".\n" ;
		return -1 ;
	}

    return 0 ;
}
 

