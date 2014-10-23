
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include "genfile/FileUtils.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "statfile/DelimitedStatSource.hpp"
#include "statfile/read_values.hpp"
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "qctool_version_autogenerated.hpp"

namespace globals {
	std::string const program_name = "checkvcf" ;
	std::string const program_version = qctool_revision ;
}

struct FixVcfOptionProcessor: public appcontext::CmdLineOptionProcessor
{
public:
	std::string get_program_name() const { return globals::program_name ; }

	void declare_options( appcontext::OptionProcessor& options ) {
		// Meta-options
		options.set_help_option( "-help" ) ;

		// File options
		options.declare_group( "Input file options" ) ;
	    options[ "-manifest" ]
	        .set_description( "Path of manifest file to compare VCF with." )
			.set_takes_values( 1 )
			.set_is_required() ;
	    options[ "-o" ]
	        .set_description( "Path of output file." )
			.set_takes_values( 1 )
		;
		options[ "-ignore-alleles" ]
			.set_description( "Do not check correctness of alleles" ) ;
		options.declare_group( "Other options" ) ;
		options [ "-log" ]
			.set_description( "Specify that " + globals::program_name + " should write a log file to the given file." )
			.set_takes_single_value() ;
		options[ "-continue-on-error" ]
			.set_description( "Specify that " + globals::program_name + " should continue if errors are found (lines with errors are skipped from output)." )
		;
	}
} ;

struct FixVCFApplication: public appcontext::ApplicationContext
{
public:
	FixVCFApplication( int argc, char** argv ):
		appcontext::ApplicationContext(
			globals::program_name,
			globals::program_version,
			std::auto_ptr< appcontext::OptionProcessor >( new FixVcfOptionProcessor ),
			argc,
			argv,
			"-log"
		)
	{
		if( options().check( "-o" ) && boost::filesystem::exists( options().get< std::string >( "-o" ) ) ) {
			get_ui_context().logger() << "Output file \"" <<  options().get< std::string >( "-o" ) << "\" exists, I will not overwrite.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		loadManifest( options().get< std::string >( "-manifest" ) ) ;
		summarise() ;
		process() ;
	}

private:
	
	typedef boost::unordered_map< std::string, genfile::SNPIdentifyingData > VariantMap ;
	VariantMap m_map ;
	
	void loadManifest( std::string const& filename ) {
		std::string const& wantedColumns = "IlmnID|Name|Chr|MapInfo|SNP" ;
		statfile::DelimitedStatSource::UniquePtr source(
			new statfile::DelimitedStatSource(
				filename,
				",",
				std::string( "[Assay]" ), std::string( "[Controls]" )
			)
		) ;
		std::map< std::size_t, std::size_t > const wantedColumnMap = statfile::impl::get_indices_of_columns(
			source->column_names(),
			wantedColumns
		) ;
		std::string IlmnID ;
		std::string Name ;
		std::string Chr ;
		genfile::Position position ;
		std::string alleleString ;
		
		appcontext::UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading manifest variants" ) ;
		progress_context( 0, source->number_of_rows() ) ;
		
		while(
			statfile::read_values(
				*source,
				wantedColumnMap,
				boost::tie( IlmnID, Name, Chr, position, alleleString )
			)
		) {
			// parse alleles
			if( alleleString.size() != 5 || alleleString[0] != '[' || alleleString[4] != ']' || alleleString[2] != '/' ) {
				throw genfile::MalformedInputError(
					"FixVCFApplication::loadManifest()",
					"Alleles for SNP " + IlmnID + " (\"" + alleleString + "\") appear malformed.",
					source->number_of_rows_read()
				) ;
			}
			std::string const allele1 = alleleString.substr( 1, 1 ) ;
			std::string const allele2 = alleleString.substr( 3, 1 ) ;
			
			genfile::SNPIdentifyingData snp(
				IlmnID, Name, genfile::GenomePosition( Chr, position ), allele1, allele2
			) ;
			
			VariantMap::const_iterator where = m_map.find( Name ) ;
			if( where != m_map.end() ) {
				throw genfile::DuplicateEntryError(
					"FixVCFApplication::loadManifest()",
					filename,
					Name,
					"Name",
					"SNP with name \"" + Name + "\" appears twice in manifest." ) ;
			}
			m_map[ Name ] = snp ;
			(*source) >> statfile::ignore_all() ;
			progress_context( source->number_of_rows_read(), source->number_of_rows() ) ;
		}
	}
	
	void summarise() const {
		get_ui_context().logger() << std::string( 72, '=' ) << "\n\n" ;

		get_ui_context().logger() << "Loaded manifest (\"" + options().get< std::string >( "-manifest" ) + "\") with the following properties:\n" ;
		get_ui_context().logger() << " -- Number of variants: " << m_map.size() << ".\n" ;
		get_ui_context().logger() << " -- First few variants: " << m_map.size() << ".\n" ;
		{
			VariantMap::const_iterator i = m_map.begin() ;
			VariantMap::const_iterator end_i = m_map.end() ;
			for( std::size_t count = 0; count < std::min( 10ul, m_map.size() ) && i != end_i; ++i, ++count ) {
				get_ui_context().logger() << " -- \"" << i->first << "\": " << i->second << ".\n" ;
			}
			get_ui_context().logger() << "\n" ;
		}
		
		get_ui_context().logger() << std::string( 72, '=' ) << "\n\n" ;
	}

	void process() {
		try {
			unsafe_process() ;
		}
		catch( genfile::InputError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::FileNotFoundError const& e ) {
			get_ui_context().logger() << "\nError: No file matching \"" << e.filespec() << "\" could be found.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}

	void unsafe_process() {
		using genfile::string_utils::to_string ;
		
		std::string line ;
		std::istream* inStream = &(std::cin) ;
		std::auto_ptr< std::ostream > outStream ;
		
		if( options().check( "-o" )) {
			outStream = genfile::open_text_file_for_output( options().get< std::string >( "-o" ) ) ;
		}
		
		appcontext::UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Processing vcf file" ) ;
		progress_context( 0, boost::optional< std::size_t >() ) ;
		
		// Ignore metadata
		std::size_t lineCount = 0 ;
		while( std::getline( *inStream, line ) && line.size() > 1 && line[0] == '#' && line[1] == '#' ) {
			if( outStream.get() ) {
				outStream->write( line.data(), line.size() ).put( '\n' ) ;
			}
			++lineCount ;
		}
		
		if( line.substr( 0, 6 ) != "#CHROM" ) {
			throw genfile::MalformedInputError( "FixVCFApplication::process()", "vcf has malformed header line", lineCount ) ;
		}
		if( outStream.get() ) {
			outStream->write( line.data(), line.size() ).put( '\n' ) ;
		}
		std::vector< std::size_t > fieldPos( 6, 0 ) ;
		bool continueOnError = options().check( "-continue-on-error" ) ;
		bool compareAlleles = !options().check( "-ignore-alleles" ) ;

		std::size_t totalVariantCount = 0 ;
		std::size_t goodVariantCount = 0 ;

		for( ; std::getline( *inStream, line ); ++lineCount, ++totalVariantCount ) {
			try {
				// otherwise look up SNP by name
				// id is third column
				for( std::size_t i = 1; i < fieldPos.size(); ++i ) {
					fieldPos[i] = line.find( '\t', fieldPos[i-1] ) ; 
					if( fieldPos[i] == std::string::npos ) {
						throw genfile::MalformedInputError(
							"FixVCFApplication::process()",
							"vcf line has fewer than " + genfile::string_utils::to_string( i+1) + " tab characters",
							lineCount
						) ;
					}
					++fieldPos[i] ; // skip the delimiter
				}

				std::string rsid = line.substr( fieldPos[2], fieldPos[3] - fieldPos[2] - 1 ) ;
				if( rsid.size() > 0 ) {
					// Just take the first ID.
					rsid = genfile::string_utils::slice( rsid ).split( "," )[0] ;
				}
				genfile::Chromosome chromosome( line.substr( fieldPos[0], fieldPos[1] - fieldPos[0] - 1 ) ) ;
				genfile::Position position = genfile::string_utils::to_repr< genfile::Position >( line.substr( fieldPos[1], fieldPos[2] - fieldPos[1] - 1 ) ) ;

				VariantMap::const_iterator where = m_map.find( rsid ) ;
				if( where == m_map.end() ) {
					throw genfile::MalformedInputError(
						"FixVCFApplication::process()",
						"vcf contains id (" + rsid + ") that is not in the variant map.",
						lineCount
					) ;
				}
				if( chromosome != where->second.get_position().chromosome() ) {
					throw genfile::MalformedInputError(
						"FixVCFApplication::process()",
						"Chromosome for SNP " + rsid + " (" + to_string( chromosome ) + ") does not match manifest chromosome (" + to_string( where->second.get_position().chromosome() ) + ")",
						lineCount
					) ;
				}
				if( position != where->second.get_position().position() ) {
					throw genfile::MalformedInputError(
						"FixVCFApplication::process()",
						"Position for SNP " + rsid + " (" + to_string( position ) + ") does not match manifest chromosome (" + to_string( where->second.get_position().position() ) + ")",
						lineCount
					) ;
				}
				if( compareAlleles ) {
					if( line.compare( fieldPos[3], fieldPos[4] - fieldPos[3] - 1, where->second.get_first_allele() ) != 0 ) {
						throw genfile::MalformedInputError(
							"FixVCFApplication::process()",
							"REF allele for SNP " + rsid + " (\"" + line.substr( fieldPos[3], fieldPos[4] - fieldPos[3] - 1 )
								+ "\") does not match manifest allele (\"" + where->second.get_first_allele() + "\")",
							lineCount
						) ;
					}
					if( line.compare( fieldPos[4], fieldPos[5] - fieldPos[4] - 1, where->second.get_second_allele() ) != 0 ) {
						throw genfile::MalformedInputError(
							"FixVCFApplication::process()",
							"REF allele for SNP " + rsid + " (\"" + line.substr( fieldPos[4], fieldPos[5] - fieldPos[4] - 1 )
								+ "\") does not match manifest allele (\"" + where->second.get_second_allele() + "\")",
							lineCount
						) ;
					}
				}
				rsid += "," + where->second.get_SNPID() ;

				if( outStream.get() ) {
					outStream->write( line.data(), fieldPos[2] ) ;
					outStream->write( rsid.data(), rsid.size() ) ;
					outStream->put( '\t' ) ;
					outStream->write( line.data() + fieldPos[3], line.size() - fieldPos[3] ) ;
					outStream->put( '\n' ) ;
				}

				++goodVariantCount ;
			}
			catch( genfile::MalformedInputError const& e ) {
				get_ui_context().logger() << "!! Error (" << e.what() << "): " << e.format_message() << ".\n" ;
				if( !continueOnError ) {
					throw ;
				}
			}
			
			progress_context( lineCount, boost::optional< std::size_t >() ) ;
		}
		progress_context( lineCount, lineCount ) ;

		get_ui_context().logger() << "Visited " << totalVariantCount << " variants of which " << goodVariantCount << " matched SNPs in the manifest.\n" ;
	}	
} ;

int main( int argc, char** argv ) {
    try {
		FixVCFApplication app( argc, argv ) ;
    }
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}
