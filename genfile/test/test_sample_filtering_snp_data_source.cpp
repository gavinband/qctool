
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <fstream>
#include "test_case.hpp"

#if HAVE_BOOST_FILESYSTEM
	#include <boost/filesystem/operations.hpp>
#endif

#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceRack.hpp"
#include "genfile/GenFileSNPDataSource.hpp"
#include "genfile/SampleFilteringSNPDataSource.hpp"
#include "stdint.h"

AUTO_TEST_SUITE( SampleFilteringSNPDataSourceTest )

// The following section contains a simple snp block writer.
namespace data {
	namespace {
		// this data has 1504 samples per row.
		unsigned int const number_of_samples = 1504 ;
	}
}

namespace {
	// The following section defines the needed objects for use with the bgen.hpp implementation.
	template< typename T >
	struct Setter
	{
		Setter( T& field ): m_field( field ) {} ;
		template< typename T2 >
		void operator()( T2 const& value ) { m_field = T( value ) ; }
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
			m_probabilities[i].AA = aa ;  
			m_probabilities[i].AB = ab ;  
			m_probabilities[i].BB = bb ;  
		}

	private:
		std::vector<probabilities_t>& m_probabilities ;
	} ;

	struct SnpData {
	
		SnpData(): probabilities( data::number_of_samples ) {} ;
	
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


	void create_file( std::string const& data, std::string const& filename ) {
		// set up our original file.
		std::ofstream file( filename.c_str() ) ;
		file << data ;
	}

	std::vector< SnpData > read_snp_data( genfile::SNPDataSource& snp_data_source ) {
		std::vector< SnpData > result ;
		SnpData snp_data ;
	
		while( snp_data_source.read_snp(
			make_setter( snp_data.number_of_samples ),
			make_setter( snp_data.SNPID ),
			make_setter( snp_data.RSID ),
			make_setter( snp_data.chromosome ),
			make_setter( snp_data.SNP_position ),
			make_setter( snp_data.allele1 ), 
			make_setter( snp_data.allele2 ),
			ProbabilitySetter( snp_data.probabilities )
		)) {
			result.push_back( snp_data ) ;
		}
		return result ;
	}

	std::vector< std::string > construct_data() {
		std::vector< std::string > data ;
		data.push_back( "" ) ;
		data.push_back( "SNP1 RS1 1000 a g 1 2 3 \n" ) ;
		data.push_back( "SNP1 RS1 1000 a g 1 2 3 4 5 6\n" ) ;
		data.push_back( "SNP1 RS1 1000 a g 1 2 3 4 5 6 7 8 9 \n" ) ;
		data.push_back( "SNP1 RS1 1000 a g 1 2 3 4 5 6 7 8 9 10 11 12\n" ) ;
		data.push_back( "SNP1 RS1 1000 a g 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15\n" ) ;
		return data ;
	}
}

AUTO_TEST_CASE( test_sample_filtering_snp_data_source ) {
	std::vector< std::string > data = construct_data() ;
	std::string const filename = tmpnam(0) + std::string( ".gen" ) ;
	for( std::size_t data_i = 1; data_i < data.size(); ++data_i ) {
		std::cerr << "====== Looking at data " << data_i << " ======\n" ;
		create_file( data[data_i], filename ) ;
		
		std::size_t number_of_samples = data_i ;
		std::set< std::set< std::size_t > > subsets ;
		// insert empty set
		subsets.insert( std::set< std::size_t >() ) ;
		// recursively construct all other sets.
		while( subsets.size() < std::pow( 2.0, double( number_of_samples ) )) {
			std::set< std::set< std::size_t > >::const_iterator
				i = subsets.begin(),
				end_i = subsets.end() ;
			for( ; i != end_i; ++i ) {
				for( std::size_t elt = 0; elt < number_of_samples; ++elt ) {
					std::set< std::size_t > new_set = *i ;
					new_set.insert( elt ) ;
					subsets.insert( new_set ) ;
				}
			}
		}

		std::set< std::set< std::size_t > >::const_iterator
			i = subsets.begin(),
			end_i = subsets.end() ;
		for( ; i != end_i; ++i ) {
			genfile::SNPDataSource::UniquePtr source( new genfile::GenFileSNPDataSource( filename, genfile::Chromosome() )) ;
			genfile::SampleFilteringSNPDataSource::UniquePtr filtering_source = genfile::SampleFilteringSNPDataSource::create(
				source,
				*i
			) ;
			
			std::size_t const filtered_number_of_samples = number_of_samples - i->size() ;
			TEST_ASSERT( filtering_source->number_of_samples() == filtered_number_of_samples ) ;

			std::vector< SnpData > data = read_snp_data( *filtering_source ) ;
			TEST_ASSERT( data.size() == 1 ) ;
			for( std::size_t j = 0; j < data.size(); ++j ) {

				TEST_ASSERT( data[j].number_of_samples == filtered_number_of_samples ) ;

				for( std::size_t elt = 0, filtered_elt = 0; elt < number_of_samples; ++elt ) {
					if( i->find( elt ) == i->end() ) {
						TEST_ASSERT( data[ j ].probabilities[ filtered_elt ].AA == (3*elt) + 1 ) ;
						TEST_ASSERT( data[ j ].probabilities[ filtered_elt ].AB == (3*elt) + 2 ) ;
						TEST_ASSERT( data[ j ].probabilities[ filtered_elt ].BB == (3*elt) + 3 ) ;
						++filtered_elt ;
					}
				}
			}
		}
	}
	std::cout << "==== success ====\n" ;
}

AUTO_TEST_SUITE_END()
