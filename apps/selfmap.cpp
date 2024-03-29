
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <map>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "genfile/MissingValue.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "components/SNPSummaryComponent/GenomeSequence.hpp"
#include "genfile/Error.hpp"
#include "genfile/db/Connection.hpp"
#include "genfile/db/SQLite3Connection.hpp"
#include "genfile/db/SQLite3Statement.hpp"
#include "genfile/db/Error.hpp"

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
				.set_takes_values_until_next_option() ;
			;
		}
		{
			options.declare_group( "Analysis options" ) ;
			options[ "-range" ]
				.set_description( "Specify one or more ranges of SNPs to operate on. "
					"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
					"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
					"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
				.set_takes_values_until_next_option()
			;
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
			options[ "-table-prefix" ]
				.set_description( "Specify a prefix for the output table names." )
				.set_takes_single_value()
				.set_default_value( "" ) ;
		}
		{
			options.declare_group( "Options affecting output" ) ;
			options[ "-o" ]
				.set_description( "Specify the path to the output file." )
				.set_takes_single_value()
				.set_is_required() ;
			options[ "-include-diagonal" ]
				.set_description( "Specify whether selfmap should include uniquely-mapping kmers in the output. "
					"To reduce the size of the output, these are not output by default." )
				;
			options[ "-log" ]
				.set_description( "Specify the path of a log file; all screen output will be copied to the file." )
				.set_takes_single_value() ;
			options[ "-omit-sequence" ]
				.set_description(
					"Do not output sequence of kmers, just the counts. "
					"(This reduces the output file size but means kmer sequence will need "
					"to be looked up again if needed for downstream analysis.)")
			;
				
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
		m_sequences(
			GenomeSequence::create(
				options().get_values< std::string >( "-sequence" ),
				get_ui_context().get_progress_context( "Reading sequence" )
			)
		),
		m_kmer_size( options().get< std::size_t >( "-kmer-size" ) ),
		m_number_of_distinct_kmers(0)
	{
	}
	
	std::vector< genfile::GenomePositionRange > parse_ranges( std::vector< std::string > const& spec ) {
		std::vector< genfile::GenomePositionRange > result ;
		result.reserve( spec.size() ) ;
		for( std::size_t i = 0; i < spec.size(); ++i ) {
			result.push_back( genfile::GenomePositionRange::parse( spec[i] ) ) ;
		}
		return result ;
	}
    
public:
	void pre_summarise() {
		get_ui_context().logger() << "SelfMapApplication: running with the following parameters:\n"
			<< "Sequence: " << m_sequences->get_summary( "  " ) << "\n"
			<< "Ranges: " ;
		for( std::size_t i = 0; i < m_ranges.size(); ++i ) {
			get_ui_context().logger() << ( i > 0 ? " " : "" ) << m_ranges[i] ;
		}
		get_ui_context().logger() << "\n" << "kmer size: " << m_kmer_size << "bp.\n" ;
	}

	void run() {
		try {
			unsafe_run() ;
		}
		catch( genfile::InputError const& e ) {
			get_ui_context().logger() << "!! (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::db::StatementPreparationError const& e ) {
			get_ui_context().logger() << "!! (" << e.what() << "): " << e.description() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::db::StatementStepError const& e ) {
			get_ui_context().logger() << "!! (" << e.what() << "): " << e.description() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}

	void post_summarise() {
	}

private:
	GenomeSequence::UniquePtr m_sequences ;
	std::size_t m_kmer_size ;
	std::vector< genfile::GenomePositionRange > m_ranges ;
	std::size_t m_number_of_distinct_kmers ;

	typedef std::pair< genfile::GenomePosition, char > KmerPosition ;
	typedef std::vector< KmerPosition > KmerPositions ;
	typedef std::pair< int, KmerPositions > KmerCountAndPositions ;
	//typedef std::map< std::string, KmerPositions > KmerMap ;
	typedef boost::unordered_map< std::string, KmerCountAndPositions > KmerMap ;
	//typedef boost::unordered_map< std::string, std::vector< std::pair< genfile::Position, char > > > KmerMap ;
	typedef boost::bimap< std::string, std::pair< genfile::GenomePosition, char > > KmerBimap ;
	
private:

	void unsafe_run() {
		if( options().check( "-range" )) {
			m_ranges = parse_ranges( options().get_values< std::string >( "-range" ) ) ;
		} else {
			m_ranges = m_sequences->get_ranges() ;
		}

		KmerMap kmerMap ;
        for( std::size_t i = 0; i < m_ranges.size(); ++i ) {
            process_range( m_ranges[i], &kmerMap ) ;
        }
		write_output( kmerMap, open_storage() ) ;
    }

    void process_range(
        genfile::GenomePositionRange const& genomeRange,
        KmerMap* result
    ) {
        using genfile::string_utils::to_string ;
        
        genfile::GenomePosition const& start = genomeRange.start() ;
        genfile::GenomePosition const& end = genomeRange.end() ;
		if( end.position() <= start.position() + m_kmer_size ) {
			get_ui_context().logger() << "!! No kmers of length " << m_kmer_size << " in the range (of length " << ( end.position() - start.position() + 1 ) << ")\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		std::string kmer ;
		std::string kmer_reverse_complement ;

		GenomeSequence::PhysicalSequenceRange range = m_sequences->get_sequence( start.chromosome(), start.position(), end.position() ) ;
		genfile::Position rangePosition = range.first.start().position() ;
		GenomeSequence::ConstSequenceIterator kmer_start = range.second.first ;
		GenomeSequence::ConstSequenceIterator kmer_end = kmer_start + m_kmer_size ;

		{
			appcontext::UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Processing kmers (" + to_string( genomeRange ) + ")" ) ;

			for( ; kmer_end != range.second.second; ++kmer_start, ++kmer_end, ++rangePosition ) {
				progress_context.notify_progress( rangePosition, range.first.end().position() - m_kmer_size ) ;
	#if DEBUG_SELFMAP
				std::cerr << "SelfMapApplication:run(): looking at kmer: " << kmer << "...\n" ;
	#endif
				kmer.assign( kmer_start, kmer_end ) ;
				genfile::string_utils::to_upper( &kmer ) ;
				kmer_reverse_complement.assign( kmer_start, kmer_end ) ;
				impl::reverse_complement( &kmer ) ;
				char orientation = '+' ;
				if( kmer.find( 'N' ) == std::string::npos ) {
					KmerPosition kmer_position(
						genfile::GenomePosition( start.chromosome(), rangePosition ),
						orientation
					) ;
					if( kmer_reverse_complement < kmer ) {
						orientation = '-' ;
						kmer.swap( kmer_reverse_complement ) ;
					}
					KmerMap::iterator where = result->find( kmer ) ;
					if( where == result->end() ) {
						result->insert(
							std::make_pair(
								kmer,
								KmerCountAndPositions(
									m_number_of_distinct_kmers++,
									KmerPositions( 1, kmer_position )
								)
							)
						) ;
					} else {
						where->second.second.push_back( kmer_position ) ;
					}
				}
			}
			progress_context.notify_progress( rangePosition, range.first.end().position() - m_kmer_size ) ;
		}
			
#if DEBUG_SELFMAP
			std::cerr << "SelfMapApplication:run(): kmerMap.size() = " << result->size() << ".\n" ;
#endif
	}

	struct Storage {
	public: 
		typedef std::unique_ptr< Storage > UniquePtr ;
		typedef genfile::db::Connection::RowId KmerId ;

		static UniquePtr create( std::string const& filename, std::string const& prefix ) {
			return UniquePtr( new Storage( filename, prefix )) ;
		}

		Storage( std::string const& filename, std::string const& table_prefix ):
			m_filename( filename ),
			m_prefix( table_prefix ),
			m_connection( genfile::db::Connection::create( filename ))
		{
			setup_storage() ;
		}
		
	public:
		
		void store_kmer(
			std::string const& kmer,
			KmerId const& id,
			int count
		) {
			m_insert_kmer_stmt
				->bind( 1, kmer )
				.bind( 2, id )
				.bind( 3, count )
				.step() ;
			m_insert_kmer_stmt->reset() ;
		}

		void store_overlap(
			KmerId const& kmer_id,
			int64_t const size,
			genfile::GenomePosition const& position1,
			char const orientation1,
			genfile::GenomePosition const& position2,
			char const orientation2
		) {
			m_insert_overlap_stmt
				->bind( 1, kmer_id )
				.bind( 2, size )
				.bind( 3, std::string( position1.chromosome() ))
				.bind( 4, position1.position() )
				.bind( 5, std::string( 1, orientation1 ) )
				.bind( 6, std::string( position2.chromosome() ))
				.bind( 7, position2.position() )
				.bind( 8, std::string( 1, orientation2 ) )
				.step() ;
			m_insert_overlap_stmt->reset() ;
		}
		
		void finalise() const {
			// nothing to do at present.
		}
		
	private:
		std::string const m_filename ;
		std::string const m_prefix ;
		genfile::db::Connection::UniquePtr m_connection ;
		genfile::db::Connection::StatementPtr m_insert_kmer_stmt ;
		genfile::db::Connection::StatementPtr m_insert_overlap_stmt ;

	private:
		void setup_storage() {
			// Performance options
			// Assume this program will write once, so we don't need
			// extensive ACID options.
			m_connection->run_statement( "PRAGMA locking_mode = EXCLUSIVE ;" ) ;
			m_connection->run_statement( "PRAGMA journal_mode = MEMORY ;" ) ;
			m_connection->run_statement( "PRAGMA synchronous = OFF ;" ) ;
			
			m_connection->run_statement(
				"CREATE TABLE `" + m_prefix + "Kmer` ( sequence TEXT NOT NULL PRIMARY KEY, id INT NOT NULL, `count` INT NOT NULL ) WITHOUT ROWID ;"
			) ;
			
			m_connection->run_statement(
				"CREATE TABLE IF NOT EXISTS `" + m_prefix + "Overlap` ( "
					"kmer_id INT NOT NULL REFERENCES `" + m_prefix + "Kmer`(id), "
					"k INT NOT NULL, "
					"chromosome1 TEXT NOT NULL, position1 TEXT NOT NULL, orientation1 TEXT NOT NULL, "
					"chromosome2 TEXT NOT NULL, position2 TEXT NOT NULL, orientation2 TEXT NOT NULL, "
					"PRIMARY KEY (chromosome1, position1, orientation1, chromosome2, position2, orientation2 ) "
				") WITHOUT ROWID ;"
			) ;
			
			m_insert_kmer_stmt = m_connection->get_statement(
				"INSERT INTO `" + m_prefix + "Kmer` (sequence, id, count) VALUES( ?, ?, ? ) ;"
			) ;

			m_insert_overlap_stmt = m_connection->get_statement(
				"INSERT INTO `" + m_prefix + "Overlap` VALUES( ?, ?, ?, ?, ?, ?, ?, ? ) ;"
			) ;
		}
	} ;

	Storage::UniquePtr open_storage() const {
		return Storage::create(
			options().get< std::string >( "-o" ),
			options().get< std::string >( "-table-prefix" )
		) ;
	}

	void write_output( KmerMap const& kmerMap, Storage::UniquePtr storage ) {

		bool const include_diagonal = options().check( "-include-diagonal" ) ;
		bool const output_sequence = !options().check( "-omit-sequence" ) ;
		using genfile::string_utils::to_string ;
		//genfile::GenomePosition position ( m_range.start().chromosome(), 0 ) ;

		{
			appcontext::UIContext::ProgressContext progress_context
				= get_ui_context().get_progress_context( "Storing kmers" ) ;
			KmerMap::const_iterator kmer_i = kmerMap.begin() ;
			KmerMap::const_iterator end_i = kmerMap.end() ;
			std::size_t count = 0 ;
			for( ; kmer_i != end_i; ++kmer_i, ++count ) {
				std::string const& kmer = kmer_i->first ;
				KmerCountAndPositions const& count_and_positions = kmer_i->second ;
				KmerPositions const& positions = count_and_positions.second ;

				Storage::KmerId const kmer_id = count_and_positions.first ;
				storage->store_kmer(
					output_sequence ? kmer : ".", kmer_id, positions.size()
				) ;
				progress_context.notify_progress( count+1, kmerMap.size() ) ;
			}
		}

		{	
			appcontext::UIContext::ProgressContext progress_context
				= get_ui_context().get_progress_context( "Storing overlaps" ) ;
			KmerMap::const_iterator kmer_i = kmerMap.begin() ;
			KmerMap::const_iterator end_i = kmerMap.end() ;
			std::size_t count = 0 ;
			for( ; kmer_i != end_i; ++kmer_i, ++count ) {
				std::string const& kmer = kmer_i->first ;
				KmerCountAndPositions const& count_and_positions = kmer_i->second ;
				KmerPositions const& positions = count_and_positions.second ;
				Storage::KmerId const kmer_id = count_and_positions.first ;

				for( std::size_t j = 0; j < positions.size(); ++j ) {
					genfile::GenomePosition const& position1 = positions[j].first ;
					char const orientation1 = positions[j].second ;
					for( std::size_t k = ( j + ( include_diagonal ? 0 : 1 ) ); k < positions.size(); ++k ) {
						genfile::GenomePosition const& position2 = positions[k].first ;
						char const orientation2 = positions[k].second ;
						storage->store_overlap(
							kmer_id,
							kmer.size(),
							position1,
							orientation1,
							position2,
							orientation2
						) ;
					}
				}
				progress_context.notify_progress( count+1, kmerMap.size() ) ;
			}
		}
		storage->finalise() ;
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
