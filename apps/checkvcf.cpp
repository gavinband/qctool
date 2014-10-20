
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/FileUtils.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPIdentifyingData.hpp"
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
		
		options.declare_group( "Other options" ) ;
		options [ "-log" ]
			.set_description( "Specify that " + globals::program_name + " should write a log file to the given file." )
			.set_takes_single_value() ;
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
				filename, ",",
				std::string( "[Assay]" ), std::string( "[Control]" )
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
		std::string line ;
		std::istream* inStream = &(std::cin) ;
		std::ostream* outStream = &(std::cout) ;
		
		// Ignore metadata
		std::size_t lineCount = 0 ;
		while( std::getline( *inStream, line ) && line.size() > 1 && line[0] == '#' && line[1] == '#' ) {
			outStream->write( line.data(), line.size() ).put( '\n' ) ;
			++lineCount ;
		}
		
		if( line.substr( 0, 6 ) != "#CHROM" ) {
			throw genfile::MalformedInputError( "FixVCFApplication::process()", "vcf has malformed header line", lineCount ) ;
		}
		outStream->write( line.data(), line.size() ).put( '\n' ) ;
		
		for( ; std::getline( *inStream, line ); ++lineCount ) {
			// otherwise look up SNP by name
			// id is third column
			std::size_t pos1 = line.find( "\t", 0 ) ;
			if( pos1 == std::string::npos ) {
				throw genfile::MalformedInputError( "FixVCFApplication::process()", "vcf line has no tab characters", lineCount ) ;
			}
			std::size_t pos2 = line.find( "\t" ) ;
			if( pos2 == std::string::npos ) {
				throw genfile::MalformedInputError( "FixVCFApplication::process()", "vcf line has only one tab character", lineCount ) ;
			}
			std::size_t pos3 = line.find( "\t" ) ;
			if( pos3 == std::string::npos ) {
				throw genfile::MalformedInputError( "FixVCFApplication::process()", "vcf line has only two tab characters", lineCount ) ;
			}
			std::string rsid = line.substr( pos2 + 1, pos3 - pos2 - 1 ) ;
			VariantMap::const_iterator where = m_map.find( rsid ) ;
			if( where == m_map.end() ) {
				throw genfile::MalformedInputError(
					"FixVCFApplication::process()",
					"vcf contains id (" + rsid + ") that is not in the variant map.",
					lineCount
				) ;
			}
			rsid += "," + where->second.get_SNPID() ;

			outStream->write( line.data(), pos2 ) ;
			outStream->write( rsid.data(), rsid.size() ) ;
			outStream->write( line.data() + pos3, line.size() - pos3 ) ;
			outStream->put( '\n' ) ;
		}
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
