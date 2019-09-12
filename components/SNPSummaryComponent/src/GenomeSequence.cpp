
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
#include "components/SNPSummaryComponent/GenomeSequence.hpp"

GenomeSequence::UniquePtr GenomeSequence::create( std::vector< std::string > const& fasta_filenames, ProgressCallback callback ) {
	return GenomeSequence::UniquePtr(
		new GenomeSequence( fasta_filenames, callback )
	) ;
}

GenomeSequence::GenomeSequence( std::string const& fasta_filename, ProgressCallback callback ):
	m_fasta_filenames( 1, fasta_filename )
{
	load_sequence( m_fasta_filenames, callback ) ;
}

GenomeSequence::GenomeSequence( std::vector< std::string > const& fasta_filenames, ProgressCallback callback ):
	m_fasta_filenames( fasta_filenames )
{
	load_sequence( fasta_filenames, callback ) ;
}

std::string const GenomeSequence::get_spec() const {
	return "GenomeSequence(" + genfile::string_utils::join( m_fasta_filenames, "\",\"" ) + ")" ;
}

void GenomeSequence::load_sequence(
	std::vector< std::string > const& filenames,
	ProgressCallback callback
) {
	for( std::size_t i = 0; i < filenames.size(); ++i ) {
		std::vector< std::string > const elts = genfile::string_utils::split( filenames[i], "=" ) ;
		if( elts.size() != 2 ) {
			throw genfile::BadArgumentError(
				"GenomeSequence::load_sequence()",
				"filenames[" + genfile::string_utils::to_string(i) + "]=\"" + filenames[i] + "\"",
				"Filespec should be in the form: name=<path to fasta file>"
			) ;
		}
		load_sequence( elts[0], elts[1] ) ;
		if( callback ) {
			callback( i + 1, filenames.size() ) ;
		}
	}
}

void GenomeSequence::load_sequence(
	std::string const& name,
	std::string const& filename
) {
	using genfile::string_utils::to_string ;

	// read header
	std::auto_ptr< std::istream > stream( genfile::open_text_file_for_input( filename )) ;
	if( !*stream ) {
		throw genfile::MalformedInputError( filename, 0 ) ;
	}

	// Ok, now read the data
	std::size_t const BUFFER_SIZE = 1024*1024 ;
	std::vector< char > buffer( BUFFER_SIZE ) ;
	std::size_t line_count = 0 ;
	std::string sequenceName = "" ;
	enum State { eHeader = 0, eSequence = 1 } ;
	State state = eHeader ;
	ChromosomeSequence sequence ;
	// State machine
	// We are either reading a header line (which may be split into seperate reads
	// because of the buffered input or a sequence line.
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
			 if( state == eHeader ) {
				 if( sequenceName.size() == 0 ) {
					if( *p != '>' ) {
						throw genfile::MalformedInputError(
							filename,
							"File does not appear to be a FASTA file (sequence does not start with a \">\" character)",
							0
						) ;
					}
				}
				sequenceName.append( p+1, segment_end ) ;
				if( segment_end != p_end ) {
					// eol, hence we have the whole sequence name
					sequence.clear() ;
					line_count = 0 ;
					state = eSequence ;
					//std::cerr << "GenomeSequence::load_sequence(): read header: \"" << sequenceName << "\".\n" ;
				}
			} else if( state == eSequence ) {
				 sequence.insert( sequence.end(), p, segment_end ) ;
				 ++line_count ;
				 if( (segment_end == (p_end-1)) || ((segment_end < (p_end-1)) && (*(segment_end+1) == '>')) ) {
 					//std::cerr << "GenomeSequence::load_sequence(): read sequence: \"" << std::string( sequence.begin(), sequence.end() ) << "\".\n" ;
					// end of sequence.
					// Store this in our map.
					std::size_t colon = sequenceName.find( ':' ) ;
					std::size_t dash = sequenceName.find( '-' ) ;
					uint32_t sequenceStart = 1 ;
					if( colon != -1 && dash != -1 && dash > colon ) {
						try {
							sequenceStart = genfile::string_utils::to_repr< uint32_t >( sequenceName.substr( colon + 1, dash - colon - 1 )) ;
						}
						catch( genfile::string_utils::StringConversionError const& e ) {
							// ignore
						}
					}
					if( m_data.find( sequenceName ) != m_data.end() ) {
						throw genfile::DuplicateKeyError( filename, "sequence name=\"" + sequenceName + "\"" ) ;
					}
					m_data[ name + ":" + sequenceName ] = ChromosomeRangeAndSequence(
						std::make_pair( sequenceStart, sequenceStart + sequence.size() ),
						sequence
					) ;
					state = eHeader ;
					sequenceName.clear() ;
				 }
			}
			p = segment_end ;
		}
	} while( *stream ) ;
}

std::string GenomeSequence::get_summary( std::string const& prefix, std::size_t column_width ) const {
	using genfile::string_utils::to_string ;
	std::string result = prefix + "GenomeSequence ("
		+ ( m_organism ? m_organism.get() : "unknown organism" ) + ", "
		+ ( m_build ? m_build.get() : "unknown build" ) + ") "
		+ "for the following regions:" ;
	std::size_t count = 0 ;
	for( SequenceData::const_iterator i = m_data.begin(); i != m_data.end(); ++i, ++count ) {
		std::string const& identifier = m_identifiers.find( i->first )->second ;
		result += "\n" + prefix + " - chromosome " + to_string( i->first ) + " (" + identifier + "):"
			+ to_string( i->second.first.first ) + "-" + to_string( i->second.first.second )
			+  " (length " + to_string( i->second.first.second - i->second.first.first ) + ")" ;
	}
	return result ;
}

std::vector< genfile::GenomePositionRange > GenomeSequence::get_ranges() const {
	std::vector< genfile::GenomePositionRange > result ;
	SequenceData::const_iterator i = m_data.begin(), end_i = m_data.end() ;
	for( ; i != end_i; ++i ) {
		result.push_back(
			genfile::GenomePositionRange(
				i->first,
				i->second.first.first,
				i->second.first.second - 1
			)
		) ;
	}
	return result ;
}

bool GenomeSequence::has_chromosome( genfile::Chromosome const& chromosome ) const {
	SequenceData::const_iterator where = m_data.find( chromosome ) ;
	return where != m_data.end() ;
}

char GenomeSequence::get_base( genfile::GenomePosition const& position ) const {
	using namespace genfile::string_utils ;
	SequenceData::const_iterator where = m_data.find( position.chromosome() ) ;
	if( where == m_data.end() ) {
		throw genfile::BadArgumentError(
			"GenomeSequence::get_base()",
			"chromosome=\"" + to_string( position.chromosome() ) + "\"",
			"Chromosome is not in the stored sequence."
		) ;
	}
	genfile::Position const sequence_start = where->second.first.first ;
	genfile::Position const sequence_end = where->second.first.second ; // one-past-the-end
	if( position.position() >= sequence_end || position.position() < sequence_start ) {
		throw genfile::BadArgumentError(
			"GenomeSequence::get_sequence()",
			"position=" + to_string( position.position() ),
			"Position is not in the sequence (" + to_string( position.chromosome() ) + ":" + to_string( sequence_start ) + "-" + to_string( sequence_end ) + ")."
		) ;
	}
	return *(where->second.second.begin() + position.position() - sequence_start) ;
}

void GenomeSequence::get_sequence( genfile::Chromosome const& chromosome, genfile::Position start, genfile::Position end, std::deque<char>* result ) const {
	PhysicalSequenceRange range = get_sequence( chromosome, start, end ) ;
	result->resize( end - start, '.' ) ;

	std::copy(
		range.second.first, range.second.second,
		result->begin() + ( range.first.start().position() - start )
	) ;
}

GenomeSequence::PhysicalSequenceRange
GenomeSequence::get_sequence( genfile::Chromosome const& chromosome, genfile::Position start, genfile::Position end ) const {
	using namespace genfile::string_utils ;
	assert( end >= start ) ;
	SequenceData::const_iterator where = m_data.find( chromosome ) ;
	if( where == m_data.end() ) {
		throw genfile::BadArgumentError(
			"GenomeSequence::get_sequence()",
			"chromosome=\"" + to_string( chromosome ) + "\"",
			"Chromosome is not in the stored sequence."
		) ;
	}
	genfile::Position const sequence_start = where->second.first.first ;
	genfile::Position const sequence_end = where->second.first.second ; // one-past-the-end
	if( start > sequence_end || start < sequence_start || end < sequence_start || end > sequence_end ) {
		throw genfile::BadArgumentError(
			"GenomeSequence::get_sequence()",
			"start..end=" + to_string( start ) + ".." + to_string( end ),
			"Region does not overlap the sequence (" + to_string( chromosome ) + ":" + to_string( sequence_start ) + "-" + to_string( sequence_end ) + ")."
		) ;
	}
	genfile::Position actual_start = std::max( start, sequence_start ) ;
	genfile::Position actual_end = std::min( end, sequence_end ) ;
	return PhysicalSequenceRange(
		genfile::GenomePositionRange(
			chromosome, actual_start, actual_end
		),
		ConstSequenceRange(
			where->second.second.begin() + (start - sequence_start),
			where->second.second.begin() + (end - sequence_start)
		)
	) ;
}
