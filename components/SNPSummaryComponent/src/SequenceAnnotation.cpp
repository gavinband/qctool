
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <utility>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <boost/optional.hpp>
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/VariantEntry.hpp"
#include "components/SNPSummaryComponent/SequenceAnnotation.hpp"

SequenceAnnotation::SequenceAnnotation( std::string const& annotation_name, std::string const& fasta_filename, ProgressCallback callback ):
	m_annotation_name( annotation_name ),
	m_fasta_filename( fasta_filename ),
	m_filenames( genfile::wildcard::find_files_by_chromosome( fasta_filename ) ),
	m_flanking( std::make_pair( 0, 0 ))
{
	load_sequence( m_filenames, &m_sequence, callback ) ;
}

void SequenceAnnotation::load_sequence( std::vector< genfile::wildcard::FilenameMatch > const& files, Sequence* sequence, ProgressCallback callback ) {
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
			callback( i + 1, files.size() ) ;
		}
	}
}

void SequenceAnnotation::load_sequence( genfile::wildcard::FilenameMatch const& file, Chromosome* chromosome, ChromosomeSequence* sequence, std::pair< genfile::Position, genfile::Position >* limits ) {
	assert( sequence ) ;
	assert( sequence->empty() ) ;
	using genfile::string_utils::to_string ;

	// read header
	std::auto_ptr< std::istream > stream( genfile::open_text_file_for_input( file.filename() )) ;
	std::string line ;
	std::getline( *stream, line ) ;
	if( !*stream ) {
		throw genfile::MalformedInputError( file.filename(), 0 ) ;
	}
	using genfile::string_utils::slice ;
	std::vector< slice > elts = slice( line ).split( ":|" ) ;
	if( elts.size() != 6 ) {
		// There seem to be many different header lines around in FASTA format files.
		// We insist on having six fields in the header which are interpreted as:
		// 1: the organism
		// 2: a build identifier
		// 3: the chromosome
		// 4,5: the range of positions represented.
		// 6: arbitrary text.
		// This is as in ancestral sequence files downloaded from ftp.1000g.ensembl.ac.uk
		// For other files users will need to insert the headers as appropriate.

		throw genfile::MalformedInputError( file.filename(), 0, elts.size() ) ;
	}
	// delimiter seems to be ":" in data downloaded from ftp.1000g.ensembl.ac.uk.  However,
	// This: http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
	// says the delimiter is |
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

	if( m_organism.is_missing() ) {
		m_organism = organism ;
		m_build = build ;
	}
	else if( m_organism != genfile::VariantEntry( organism ) ) {
		throw genfile::MismatchError( "SequenceAnnotation::load_sequence()", m_fasta_filename, "organism=\"" + to_string( m_organism ) + "\"", "organism=\"" + to_string( organism ) + "\"" ) ;
	}
	else if( m_build != build ) {
		throw genfile::MismatchError( "SequenceAnnotation::load_sequence()", m_fasta_filename, "build=\"" + to_string( m_build ) + "\"", "build=\"" + to_string( build ) + "\"" ) ;
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

std::string SequenceAnnotation::get_summary( std::string const& prefix, std::size_t column_width ) const {
	using genfile::string_utils::to_string ;
	std::string result = prefix + "SequenceAnnotation: loaded " + m_annotation_name + " sequence (" + to_string( m_organism ) + ", build " + to_string( m_build ) + ") for the following regions:\n" ;
	using genfile::string_utils::to_string ;
	std::size_t count = 0 ;
	for( Sequence::const_iterator i = m_sequence.begin(); i != m_sequence.end(); ++i, ++count ) {
		if( m_sequence.size() > 5 && count == 1 ) {
			std::advance( i, m_sequence.size() - count - 2 ) ;
			count = m_sequence.size() - 2 ;
			result += "\n" + prefix + " - ." ;
			result += "\n" + prefix + " - ." ;
		}
		result += "\n" + prefix + " - chromosome " + to_string( i->first ) + ", length "
			+ to_string( 1 + i->second.first.second - i->second.first.first ) + " ("
			+ to_string( i->second.first.first ) + "-" + to_string( i->second.first.second )
			+ ")" ;
	}
	return result ;
}

void SequenceAnnotation::set_flanking( std::size_t left, std::size_t right ) {
	m_flanking.first = left ;
	m_flanking.second = right ;
}

void SequenceAnnotation::operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const&, genfile::VariantDataReader&, ResultCallback callback ) {
	using namespace genfile::string_utils ;
	Sequence::const_iterator where = m_sequence.find( snp.get_position().chromosome() ) ;
	if( where != m_sequence.end() ) {
		genfile::Position pos = snp.get_position().position() ;
		std::size_t const start = where->second.first.first ;
		std::size_t const end = where->second.first.second ;
		if( pos >= start && pos <= end ) {
			std::string allele( 1, where->second.second[ pos - start ] ) ;
			if( allele != "." ) {
				//std::cerr << snp << ": ancestral allele is " << std::string( 1, allele ) << ".\n";
				callback( m_annotation_name + "_allele", allele ) ;
				if( to_upper( allele ) != snp.get_first_allele() && to_upper( allele ) != snp.get_second_allele() ) {
					callback( m_annotation_name + "_allele_warning", m_annotation_name + " allele does not match the variant alleles; is strand alignment required?" ) ;
				}
			}
		
			if( m_flanking.first > 0 || m_flanking.second > 0 ) {
				std::size_t const left_flanking_start = pos - std::min< genfile::Position >( pos, m_flanking.first ) ;
				std::size_t const left_flanking_end = pos ;
				std::size_t const right_flanking_start = pos + std::min< genfile::Position >( end - pos, 1ul ) ;
				std::size_t const right_flanking_end = pos + std::min< genfile::Position >( end - pos, m_flanking.second ) ;
				// Figure out which of our alleles is the reference
				std::ostringstream flank ;
				std::copy(
					where->second.second.begin() + left_flanking_start - start,
					where->second.second.begin() + left_flanking_end - start,
					std::ostream_iterator< char >( flank )
				) ;
				flank << "[" << allele ;
				if( to_upper( snp.get_first_allele() ) == to_upper( allele ) ) {
					flank << "/" << snp.get_second_allele() ;
				}
				else if( to_upper( snp.get_second_allele() ) == to_upper( allele ) ) {
					flank << "/" << snp.get_first_allele() ;
				}
				else {
					flank << "/" << snp.get_first_allele() << "/" << snp.get_second_allele() ;
				}
				flank << "]" ;
				std::copy(
					where->second.second.begin() + right_flanking_start - start,
					where->second.second.begin() + right_flanking_end - start,
					std::ostream_iterator< char >( flank )
				) ;
					
				
				callback(
					m_annotation_name + "_flanking_sequence",
					flank.str()
				) ;
			}
		}
	}
}
