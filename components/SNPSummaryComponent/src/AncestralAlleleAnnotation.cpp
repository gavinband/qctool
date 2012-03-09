#include <fstream>
#include <boost/optional.hpp>
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils/slice.hpp"
#include "components/SNPSummaryComponent/AncestralAlleleAnnotation.hpp"

AncestralAlleleAnnotation::AncestralAlleleAnnotation( std::string const& fasta_filename, ProgressCallback callback ):
	m_fasta_filename( fasta_filename ),
	m_filenames( genfile::wildcard::find_files_by_chromosome( fasta_filename ) )
{
	load_sequence( m_filenames, &m_sequence, callback ) ;
}

void AncestralAlleleAnnotation::load_sequence( std::vector< genfile::wildcard::FilenameMatch > const& files, Sequence* sequence, ProgressCallback callback ) {
	assert( sequence ) ;
	for( std::size_t i = 0; i < files.size(); ++i ) {
		ChromosomeRangeAndSequence data ;
		Chromosome chromosome ;
		load_sequence( files[i], &chromosome, &data.second, &data.first ) ;
		if( m_sequence.find( chromosome ) != m_sequence.end() ) {
			throw genfile::DuplicateKeyError( m_fasta_filename, "chromosome=\"" + files[i].match() + "\"" ) ;
		}
		m_sequence[ chromosome ].first = data.first ;
		m_sequence[ chromosome ].second.swap( data.second ) ;

		if( callback ) {
			callback( i, boost::optional< std::size_t >() ) ;
		}
	}
}

void AncestralAlleleAnnotation::load_sequence( genfile::wildcard::FilenameMatch const& file, Chromosome* chromosome, ChromosomeSequence* sequence, std::pair< std::size_t, std::size_t >* limits ) {
	assert( sequence ) ;
	assert( sequence->empty() ) ;
	// read header
	std::auto_ptr< std::istream > stream( genfile::open_text_file_for_input( file.filename() )) ;
	std::string line ;
	std::getline( *stream, line ) ;
	if( !*stream ) {
		throw genfile::MalformedInputError( file.filename(), 0 ) ;
	}
	using genfile::string_utils::slice ;
	// delimiter seems to be ":" in data downloaded from ftp.1000g.ensembl.ac.uk.  However,
	// This: http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
	// says the delimiter is |
	std::vector< slice > elts = slice( line ).split( ":|" ) ;
	if( elts.size() != 6 ) {
		throw genfile::MalformedInputError( file.filename(), 0, elts.size() ) ;
	}

	std::string organism = elts[0] ;
	std::string build = elts[1] ;
	Chromosome read_chromosome( elts[2] ) ;
	std::pair< std::size_t, std::size_t > read_limits ;

	try {
		read_limits.first = genfile::string_utils::to_repr< std::size_t >( elts[3] ) ;
	}
	catch( genfile::string_utils::StringConversionError const& ) {
		throw genfile::MalformedInputError( file.filename(), 0, 3 ) ;
	}
	try {
		read_limits.second = genfile::string_utils::to_repr< std::size_t >( elts[4] ) ;
	}
	catch( genfile::string_utils::StringConversionError const& ) {
		throw genfile::MalformedInputError( file.filename(), 0, 4 ) ;
	}

	if( m_organism == "" ) {
		m_organism = organism ;
		m_build = build ;
	}
	else if( m_organism != organism ) {
		throw genfile::MismatchError( "AncestralAlleleAnnotation::load_sequence()", m_fasta_filename, "organism=\"" + m_organism + "\"", "organism=\"" + organism + "\"" ) ;
	}
	else if( m_build != build ) {
		throw genfile::MismatchError( "AncestralAlleleAnnotation::load_sequence()", m_fasta_filename, "build=\"" + m_build + "\"", "build=\"" + build + "\"" ) ;
	}

	(*limits) = read_limits ;
	(*chromosome) = read_chromosome ;
	
	// Ok, now read the data
	std::size_t const BUFFER_SIZE = 1024*1024 ;
	std::vector< char > buffer( BUFFER_SIZE ) ;
	std::size_t line_count = 0 ;
	do {
		stream->read( &buffer[0], BUFFER_SIZE ) ;
		std::size_t const count = stream->gcount() ;
		char const* p_end = &buffer[0] + count ;
		for(
			char const* p = &buffer[0] ;
			p < p_end;
			++p
		 ) {
			char const* segment_end = std::find( p, p_end, '\n' ) ;
			sequence->insert( sequence->end(), p, segment_end ) ;
			p = segment_end ;
			++line_count ;
		}
	} while( *stream ) ;
	
	if( sequence->size() != ( 1 + limits->second - limits->first ) ) {
		throw genfile::MalformedInputError( file.filename(), line_count ) ;
	}
}

std::string AncestralAlleleAnnotation::get_summary( std::string const& prefix, std::size_t column_width ) const {
	std::string result = prefix + "AncestralAlleleAnnotation: using the following files:\n" ;
	using genfile::string_utils::to_string ;
	for( Sequence::const_iterator i = m_sequence.begin(); i != m_sequence.end(); ++i ) {
		result += prefix + " - chromosome " + to_string( i->first ) + ", length "
			+ to_string( i->second.first.second - i->second.first.first ) + " ("
			+ to_string( i->second.first.first ) + "-" + to_string( i->second.first.second )
			+ ")\n" ;
	}
	return result ;
}

void AncestralAlleleAnnotation::operator()( SNPIdentifyingData const& snp, Genotypes const&, genfile::VariantDataReader&, ResultCallback callback ) {
	Sequence::const_iterator where = m_sequence.find( snp.get_position().chromosome() ) ;
	if( where != m_sequence.end() ) {
		genfile::Position pos = snp.get_position().position() ;
		if( pos >= where->second.first.first && pos <= where->second.first.second ) {
			char allele = where->second.second[ pos - where->second.first.first ] ;
			if( allele != '.' ) {
				//std::cerr << snp << ": ancestral allele is " << std::string( 1, allele ) << ".\n";
				callback( "ancestral allele", std::string( 1, allele ) ) ;
			}
		}
	}
}
