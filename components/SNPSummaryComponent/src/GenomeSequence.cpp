
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

GenomeSequence::UniquePtr GenomeSequence::create( std::string const& fasta_filename, ProgressCallback callback ) {
	return GenomeSequence::UniquePtr(
		new GenomeSequence( fasta_filename, callback )
	) ;
}

GenomeSequence::GenomeSequence( std::string const& fasta_filename, ProgressCallback callback ):
	m_fasta_filename( fasta_filename )
{
	load_sequence( genfile::wildcard::find_files_by_chromosome( fasta_filename, genfile::wildcard::eALL_CHROMOSOMES ), &m_sequence, callback ) ;
}

void GenomeSequence::load_sequence( std::vector< genfile::wildcard::FilenameMatch > const& files, SequenceData* sequence, ProgressCallback callback ) {
	assert( sequence ) ;
	for( std::size_t i = 0; i < files.size(); ++i ) {
		ChromosomeRangeAndSequence data ;
		Chromosome chromosome ;
		std::string identifier ;
		load_sequence( files[i], &chromosome, &data.second, &data.first, &identifier ) ;
		if( m_sequence.find( chromosome ) != m_sequence.end() ) {
			throw genfile::DuplicateKeyError( m_fasta_filename, "chromosome=\"" + files[i].match() + "\"" ) ;
		}
		m_sequence[ chromosome ].first = data.first ;
		m_sequence[ chromosome ].second.swap( data.second ) ;
		m_identifiers[ chromosome ] = identifier ;
		if( callback ) {
			callback( i + 1, files.size() ) ;
		}
	}
}

void GenomeSequence::load_sequence(
	genfile::wildcard::FilenameMatch const& file,
	Chromosome* chromosome,
	ChromosomeSequence* sequence,
	std::pair< genfile::Position, genfile::Position >* limits,
	std::string* identifier
) {
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
	if( line.size() == 0 || line[0] != '>' ) {
		throw genfile::MalformedInputError(
			file.filename(),
			"File does not appear to be a FASTA file (it does not start with a line \">....\")",
			0
		) ;
	}
	line = line.substr( 1, line.size() ) ;
	*identifier = line ;
	*chromosome = genfile::Chromosome( file.match() ) ;
	
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
	
	limits->first = 1 ;
	limits->second = sequence->size() + 1 ;
}

std::string GenomeSequence::get_summary( std::string const& prefix, std::size_t column_width ) const {
	using genfile::string_utils::to_string ;
	std::string result = prefix + "GenomeSequence ("
		+ ( m_organism ? m_organism.get() : "unknown organism" ) + ", "
		+ ( m_build ? m_build.get() : "unknown build" ) + ") "
		+ "for the following regions:" ;
	std::size_t count = 0 ;
	for( SequenceData::const_iterator i = m_sequence.begin(); i != m_sequence.end(); ++i, ++count ) {
		std::string const& identifier = m_identifiers.find( i->first )->second ;
		result += "\n" + prefix + " - chromosome " + to_string( i->first ) + " (" + identifier + "):"
			+ to_string( i->second.first.first ) + "-" + to_string( i->second.first.second )
			+  " (length " + to_string( i->second.first.second - i->second.first.first ) + ")" ;
	}
	return result ;
}

bool GenomeSequence::has_chromosome( genfile::Chromosome const& chromosome ) const {
	SequenceData::const_iterator where = m_sequence.find( chromosome ) ;
	return where != m_sequence.end() ;
}

char GenomeSequence::get_base( genfile::GenomePosition const& position ) const {
	using namespace genfile::string_utils ;
	SequenceData::const_iterator where = m_sequence.find( position.chromosome() ) ;
	if( where == m_sequence.end() ) {
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
	SequenceData::const_iterator where = m_sequence.find( chromosome ) ;
	if( where == m_sequence.end() ) {
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
			where->second.second.begin() + start - sequence_start,
			where->second.second.begin() + end - sequence_start
		)
	) ;
}
