
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <deque>
#include <iostream>
#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include "genfile/string_utils.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/ProgramFlow.hpp"
#include "config/version_autogenerated.hpp"
#include "config/config.hpp"
#if HAVE_MGL
	#include "mgl/mgl.h"
	#include "mgl/mgl_zb.h"
#endif

namespace globals {
	std::string const program_name = "inflation" ;
	std::string const program_version = qctool_revision ;
}

struct InflationApplication: public appcontext::ApplicationContext
{
public:
	static void declare_options( appcontext::OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;
		options[ "-analysis-name" ]
			.set_description( "Specify a name to label results from this analysis with.  (This applies to modules which store their results in a qcdb file.)" )
			.set_takes_single_value()
			.set_default_value( "inflation analysis, started " + appcontext::get_current_time_as_string() ) ;
		options [ "-qq-plot" ]
			.set_description( "Specify that " + globals::program_name + " should make a qq-plot in the file with the given name." )
			.set_takes_single_value() ;
		options [ "-log" ]
			.set_description( "Specify that " + globals::program_name + " should write a log file to the given file." )
			.set_takes_single_value() ;
	}

public:
	InflationApplication( appcontext::OptionProcessor::UniquePtr options, int argc, char** argv ):
		appcontext::ApplicationContext( globals::program_name, globals::program_version , options, argc, argv, "-log" )
	{
	}
	
	void process() {
		unsafe_process(); 
	}
private:
	
	void unsafe_process() {
		get_ui_context().logger() << application_name() << ": reading..." ;

		std::string line ;
		std::deque< float > numbers ;
		std::size_t missing_count = 0 ;
		while( std::getline( std::cin, line )) {
			if( line.empty() || line.find_first_not_of( "0123456789.E+-" ) != std::string::npos ) {
				++missing_count ;
			} else {
				numbers.push_back( genfile::string_utils::to_repr< float >( line ) ) ;
			}

			if( ( numbers.size() + missing_count ) % 1000000 == 0 ) {
				get_ui_context().logger() << ( numbers.size() + missing_count ) << ".." ;
			}
		}

		get_ui_context().logger() << ( numbers.size() / 1000000 ) << ".\n" ;
		get_ui_context().logger() << "inflation.cpp: read " << numbers.size() << " P-values from std::cin.\n" ;
		get_ui_context().logger() << "inflation.cpp: converting to chi-squared statistics...\n" ;
		using namespace boost::math ;
		chi_squared_distribution< float > chi_square( 1 ) ;
		for( std::size_t i = 0; i < numbers.size(); ++i ) {
			numbers[i] = quantile( complement( chi_square, numbers[i] )) ;
		}

		get_ui_context().logger() << "inflation.cpp: finding median...\n" ;

		std::cout << "analysis,N,missing,median,lambda,1pc_adjusted_median,1pc_adjusted_lambda\n" ;
		std::cout
			<< options().get< std::string >( "-analysis-name" ) << ","
			<< numbers.size() << ","
			<< missing_count ;
		// Compute basic lambda
		{
			std::size_t const mid = numbers.size() / 2 ;
			std::nth_element( numbers.begin(), numbers.begin() + mid, numbers.end() ) ;
			double const median = numbers[ mid ] ;
			std::cout << "," << median << "," << (median / quantile( complement( chi_square, 0.5 ))) ;
		}

		if( options().check( "-qq-plot" )) {
#if HAVE_MGL
			get_ui_context().logger() << "inflation.cpp: making qq-plot...\n" ;	
			plot( numbers, options().get< std::string >( "-qq-plot" ), options().get< std::string >( "-analysis-name" ) ) ;
#else
			get_ui_context().logger() << "!! InflationApplication::unsafe_process(): you must build with MathGL support for -qq-plot to be supported.  Skipping plots...\n" ;
#endif
		}

		// Compute lambda after removing top 99% of signal
		{
			std::size_t const top_one_percent = numbers.size() * 0.99 ;
			std::deque< float >::iterator top_one_percent_i = numbers.begin() + top_one_percent ;
			numbers.erase( top_one_percent_i, numbers.end() ) ;
		}

		{
			std::size_t const mid = numbers.size() / 2 ;
			std::nth_element( numbers.begin(), numbers.begin() + mid, numbers.end() ) ;
			double const median = numbers[ mid ] ;
			std::cout << "," << median << "," << (median / quantile( complement( chi_square, 0.5 ))) ;
		}
		std::cout << "\n" ;
	}
	
#if HAVE_MGL
	void plot( std::deque< float >& numbers, std::string const& filename, std::string const& title ) {
		std::sort( numbers.begin(), numbers.end() ) ;

		mglData
			observed( numbers.size() ),
			expected( numbers.size() )
		;
		
		using namespace boost::math ;
		{
			chi_squared_distribution< float > chi_squared( 1 ) ;
			for( std::size_t i = 0; i < numbers.size(); ++i ) {
				observed.a[i] = numbers[i] ;
				expected.a[i] = quantile( complement( chi_squared, float( numbers.size() - i ) / ( numbers.size() + 1 ) ) ) ;
			}
		}
		mglGraphZB graph( 800, 800 ) ;
		graph.Light( true ) ;
		graph.Clf( mglColor( 1.0, 1.0, 1.0 ) ) ;
		graph.Title( title.c_str(), 0, 4 ) ;
		graph.SetTickLen( 0.04 ) ;
		graph.SetRanges(
			0.0, 30,
			0.0, 30
		) ;
		graph.Axis( "x", true ) ;
		graph.Axis( "y", true ) ;
		graph.CAxis( 0.0, 3.0 ) ;

		graph.Plot( expected, observed ) ; //, observed ) ;
		graph.WritePNG( filename.c_str(), "", false ) ;
	}
#endif
} ;

int main( int argc, char** argv ) {
	try {
		appcontext::OptionProcessor::UniquePtr options( new appcontext::OptionProcessor ) ;
		InflationApplication::declare_options( *options ) ;
		InflationApplication app( options, argc, argv ) ;
		app.process() ;
	} catch ( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}

