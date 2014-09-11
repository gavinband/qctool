
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <map>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "genfile/MissingValue.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "components/SNPSummaryComponent/GenomeSequence.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"
#include "qcdb/Storage.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"
#include "qcdb/FlatFileOutputter.hpp"

// #define DEBUG_SELFMAP 1

namespace globals {
	std::string const program_name = "selfmap" ;
	std::string const program_version = "dev" ;
}

namespace {
	double const NA = std::numeric_limits< double >::quiet_NaN() ;
	namespace impl {
		char complement_allele( char a ) {
			switch( a ) {
				case 'A': a = 'T'; break ;
				case 'G': a = 'C'; break ;
				case 'C': a = 'G'; break ;
				case 'T': a = 'A'; break ;
			}
			return a ;
		}
	
		void reverse_complement( std::string* sequence ) {
			std::reverse( sequence->begin(), sequence->end() ) ;
			std::transform( sequence->begin(), sequence->end(), sequence->begin(), &impl::complement_allele ) ;
		}
	}
}

struct SelfMapOptions: public appcontext::CmdLineOptionProcessor {
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;

		{
			options.declare_group( "Input data options" ) ;
			options[ "-sequence" ]
				.set_description( "Specify the path of a file containing sequence to load." )
				.set_is_required()
				.set_takes_single_value() ;
			;
		}
		{
			options.declare_group( "Analysis options" ) ;
			options[ "-range" ]
				.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to operate on. "
					"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
					"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
					"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
				.set_takes_single_value()
				.set_is_required() ;
			options[ "-kmer-size" ]
				.set_is_required()
				.set_takes_single_value()
				.set_description( "Specify the size of kmer to compute using." )
			;
			options[ "-analysis-name" ]
				.set_description( "Specify a name to label results from this analysis with.  (This applies to modules which store their results in a qcdb file.)" )
				.set_takes_single_value()
				.set_default_value( "qctool analysis" ) ;
			options[ "-analysis-chunk" ]
				.set_description( "Specify a name denoting the current genomic region or chunk on which this is run.  This is intended for use in parallel environments." )
				.set_takes_single_value()
				.set_default_value( genfile::MissingValue() ) ;
			options[ "-table-name" ]
				.set_description( "Specify a name for the table to use when using -flat-table." )
				.set_takes_single_value() ;
		}
		{
			options.declare_group( "Options affecting output" ) ;
			options[ "-o" ]
				.set_description( "Specify the path to the output file." )
				.set_is_required()
				.set_takes_single_value() ;
			options[ "-omit-diagonal" ]
				.set_description( "Specify that selfmap should not output matches on the diagonal." )
					;
			options[ "-log" ]
				.set_description( "Specify the path of a log file; all screen output will be copied to the file." )
				.set_takes_single_value() ;
		}

		{
			options.declare_group( "Analysis options" ) ;
		}
	}
} ;

struct SelfMapApplication: public appcontext::ApplicationContext {
public:
	typedef std::map< std::string, std::vector< std::string > > GroupDefinition ;
	typedef std::map< std::string, std::vector< std::string > > ValueListSet ;

public:
	SelfMapApplication( int argc, char **argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version,
			std::auto_ptr< appcontext::OptionProcessor >( new SelfMapOptions ),
			argc,
			argv,
			"-log"
		),
		m_sequence( GenomeSequence::create( options().get< std::string >( "-sequence" ), get_ui_context().get_progress_context( "Reading sequence" ) ) ),
		m_kmer_size( options().get< std::size_t >( "-kmer-size" ) ),
		m_range( genfile::GenomePositionRange::parse( options().get< std::string >( "-range" ) ) )
	{
	}
	
public:
	void pre_summarise() {
		get_ui_context().logger() << "SelfMapApplication: running with the following parameters:\n"
			<< "Sequence: " << m_sequence->get_summary( "  " ) << "\n"
			<< "Range: " << m_range << "\n"
			<< "kmer size: " << m_kmer_size << "bp.\n" ;
	}
	void run() {
		try {
			unsafe_run() ;
		}
		catch( genfile::InputError const& e ) {
			get_ui_context().logger() << "!! (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( db::StatementPreparationError const& e ) {
			get_ui_context().logger() << "!! (" << e.what() << "): " << e.description() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}
	void post_summarise() {
	}
private:
	GenomeSequence::UniquePtr m_sequence ;
	std::size_t m_kmer_size ;
	genfile::GenomePositionRange m_range ;

	typedef std::map< std::string, std::vector< std::pair< genfile::Position, char > > > KmerMap ;
	//typedef boost::unordered_map< std::string, std::vector< std::pair< genfile::Position, char > > > KmerMap ;
	
private:

	void unsafe_run() {
		genfile::GenomePosition const start = m_range.start() ;
		genfile::GenomePosition const end = m_range.end() ;

		if( end.position() <= start.position() + m_kmer_size ) {
			get_ui_context().logger() << "!! No kmers of length " << m_kmer_size << " in the range (of length " << ( end.position() - start.position() + 1 ) << ")\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}

		KmerMap kmerMap ;
		std::string kmer ;
		char strand = '+' ;

		GenomeSequence::PhysicalSequenceRange range = m_sequence->get_sequence( m_range.start().chromosome(), m_range.start().position(), m_range.end().position() ) ;
		genfile::Position rangePosition = range.first.start().position() ;
		GenomeSequence::ConstSequenceIterator kmer_start = range.second.first ;
		GenomeSequence::ConstSequenceIterator kmer_end = kmer_start + m_kmer_size ;


		{
			appcontext::UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Processing kmers" ) ;

			for( ; kmer_end != range.second.second; ++kmer_start, ++kmer_end, ++rangePosition ) {
				progress_context.notify_progress( rangePosition, range.first.end().position() - m_kmer_size ) ;
	#if DEBUG_SELFMAP
				std::cerr << "SelfMapApplication:run(): looking at kmer: " << kmer << "...\n" ;
	#endif
				kmer.assign( kmer_start, kmer_end ) ;
				genfile::string_utils::to_upper( &kmer ) ;
				if( kmer.find( 'N' ) == std::string::npos ) {
					kmerMap[ kmer ].push_back( std::make_pair( rangePosition, '+' ) ) ;
					impl::reverse_complement( &kmer ) ;
					kmerMap[ kmer ].push_back( std::make_pair( rangePosition + m_kmer_size, '-' ) ) ;
				}
			}
			progress_context.notify_progress( rangePosition, range.first.end().position() - m_kmer_size ) ;
		}
			
#if DEBUG_SELFMAP
			std::cerr << "SelfMapApplication:run(): kmerMap.size() = " << kmerMap.size() << ".\n" ;
#endif
		write_output( kmerMap, open_storage() ) ;
	}

	qcdb::Storage::SharedPtr open_storage() const {
		qcdb::Storage::SharedPtr storage ;
		std::string const& filename = options().get< std::string >( "-o" ) ;

		if( filename.size() > 7 && filename.substr( filename.size() - 7, 7 ) == ".sqlite" ) {
			snp_summary_component::DBOutputter::SharedPtr outputter = snp_summary_component::DBOutputter::create_shared(
				options().get< std::string >( "-o" ),
				options().get< std::string >( "-analysis-name" ),
				options().get< std::string >( "-analysis-chunk" ),
				options().get_values_as_map()
			) ;
			if( options().check( "-table-name" )) {
				outputter->set_table_name( options().get< std::string >(  "-table-name" )) ;
			}
			storage = outputter ;
		} else {
			storage = qcdb::FlatFileOutputter::create_shared(
				options().get< std::string >( "-o" ),
				options().get< std::string >( "-analysis-name" ),
				options().get_values_as_map()
			) ;
		}
		return storage ;
	}

	void write_output( KmerMap const& kmerMap, qcdb::Storage::SharedPtr storage ) {
		appcontext::UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Storing results" ) ;

		bool const omit_diagonal = options().check( "-omit-diagonal" ) ;
		using genfile::string_utils::to_string ;
		std::string const stub = to_string( m_range.start().chromosome() ) + ":" ;
		genfile::GenomePosition position ( m_range.start().chromosome(), 0 ) ;
		KmerMap::const_iterator i = kmerMap.begin() ;
		KmerMap::const_iterator end_i = kmerMap.end() ;
		std::size_t count = 0 ;
		for( ; i != end_i; ++i, ++count ) {
			position.position() = i->second[0].first ;
			genfile::SNPIdentifyingData2 snp(
				stub + to_string( i->second[0].first ),
				".",
				position,
				i->first, "."
			) ;
			for( std::size_t j = 0; j < i->second.size(); ++j ) {
				snp.set_position( genfile::GenomePosition( position.chromosome(), i->second[j].first ) ) ;
				storage->store_per_variant_data(
					snp,
					"orientation",
					std::string( 1, i->second[j].second )
				) ;
				for( std::size_t k = ( j + ( omit_diagonal ? 1 : 0 ) ); k < i->second.size(); ++k ) {
					storage->store_per_variant_data(
						snp,
						"other_position",
						genfile::VariantEntry::Integer( i->second[k].first )
					) ;

					storage->store_per_variant_data(
						snp,
						"other_orientation",
						std::string( 1, i->second[k].second )
					) ;
				}
			}
			progress_context.notify_progress( count+1, kmerMap.size() ) ;
		}
	}
} ;

int main( int argc, char** argv ) {
	try {
		SelfMapApplication app( argc, argv ) ;	
		app.pre_summarise() ;
		app.run() ;
		app.post_summarise() ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
