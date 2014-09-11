
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
	m_sequence( fasta_filename, callback ),
	m_flanking( std::make_pair( 0, 0 ))
{}

std::string SequenceAnnotation::get_summary( std::string const& prefix, std::size_t column_width ) const {
	using genfile::string_utils::to_string ;
	std::string result = prefix + "SequenceAnnotation: loaded " + m_annotation_name + " sequence:\n"
		+ m_sequence.get_summary( prefix + "  ", column_width ) ;
	return result ;
}

void SequenceAnnotation::set_flanking( std::size_t left, std::size_t right ) {
	m_flanking.first = left ;
	m_flanking.second = right ;
}

void SequenceAnnotation::operator()( SNPIdentifyingData const& snp, Genotypes const& genotypes, SampleSexes const&, genfile::VariantDataReader&, ResultCallback callback ) {
	using namespace genfile::string_utils ;
	std::deque< char > flankingSequence ;
	try {
		{
			std::size_t const allele_size = snp.get_first_allele().size() ;
			m_sequence.get_sequence( snp.get_position().chromosome(), snp.get_position().position() - m_flanking.first, snp.get_position().position() + allele_size + m_flanking.second, &flankingSequence ) ;
			assert( flankingSequence.size() == m_flanking.first + m_flanking.second + allele_size ) ;
			callback( m_annotation_name + "_" + "alleleA", std::string( flankingSequence.begin() + m_flanking.first, flankingSequence.begin() + m_flanking.first + allele_size ) ) ;
			callback( m_annotation_name + "_" + "left_flanking", std::string( flankingSequence.begin(), flankingSequence.begin() + m_flanking.first ) ) ;
			callback( m_annotation_name + "_" + "alleleA_right_flanking", std::string( flankingSequence.begin() + m_flanking.first + allele_size, flankingSequence.end() ) ) ;
		}

		{
			std::size_t const allele_size = snp.get_second_allele().size() ;
			m_sequence.get_sequence( snp.get_position().chromosome(), snp.get_position().position() - m_flanking.first, snp.get_position().position() + allele_size + m_flanking.second, &flankingSequence ) ;
			assert( flankingSequence.size() == m_flanking.first + m_flanking.second + allele_size ) ;
			callback( m_annotation_name + "_" + "alleleB", std::string( flankingSequence.begin() + m_flanking.first, flankingSequence.begin() + m_flanking.first + allele_size ) ) ;
			callback( m_annotation_name + "_" + "alleleA_right_flanking", std::string( flankingSequence.begin() + m_flanking.first + allele_size, flankingSequence.end() ) ) ;
		}
	}
	catch( genfile::BadArgumentError const& e ) {
		// nowt to do.
	}
}
