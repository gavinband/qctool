#include "GenRow.hpp"
#include "GenRowSimpleFileSource.hpp"
#include "FileUtil.hpp"
#include "Whitespace.hpp"

// Read gen rows from a gen file making no attempt to chunk reads
GenRowSimpleFileSource::GenRowSimpleFileSource( INPUT_FILE_PTR stream_ptr )
: m_stream_ptr( stream_ptr )
{}

void GenRowSimpleFileSource::read_genrow( GenRow& row ) {
	assert( !check_if_empty() ) ;
	// read trailing whitespace so that stream's eof will flag.
	Whitespace whitespace ;
	(*m_stream_ptr) >> whitespace >> row >> whitespace ;
}

bool GenRowSimpleFileSource::check_if_empty() {
	return !m_stream_ptr.get() || !m_stream_ptr->good() ;
}

