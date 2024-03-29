
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "config/config.hpp"
#if HAVE_BOOST_FUNCTION
#include <boost/function.hpp>
#endif
#include <boost/optional.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/SNPDataSourceRack.hpp"
#include "genfile/wildcard.hpp"

namespace genfile {
	SNPDataSourceChain::UniquePtr SNPDataSourceChain::create(
		std::vector< wildcard::FilenameMatch > const& filenames,
		boost::optional< vcf::MetadataParser::Metadata > const& metadata,
		std::string const& filetype_hint,
		NotifyProgress notify_progress
	) {
		std::auto_ptr< SNPDataSourceChain > chain( new SNPDataSourceChain() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			if( notify_progress ) {
				notify_progress( i, filenames.size() ) ;
			}
			Chromosome chromosome ;
			if( filenames[i].match() != "" ) {
				chromosome = filenames[i].match() ;
			}
			chain->add_source( SNPDataSource::create( filenames[i].filename(), chromosome, metadata, filetype_hint )) ;
			if( notify_progress ) {
				notify_progress( i + 1, filenames.size() ) ;
			}
		}
		return chain ;
	}

	SNPDataSourceChain::UniquePtr SNPDataSourceChain::create(
		std::vector< std::vector< wildcard::FilenameMatch > > const& filenames,
		NotifyProgress notify_progress
	) {
		for( std::size_t i = 1; i < filenames.size(); ++i ) {
			assert( filenames[i].size() == filenames[0].size() ) ;
		}

		std::size_t total_number_of_files = filenames.empty() ? 0 : ( filenames[0].size() * filenames.size() ) ;
		std::size_t count = 0 ;

		std::auto_ptr< SNPDataSourceChain > chain( new SNPDataSourceChain() ) ;
		if( filenames.size() > 0 ) {
			for( std::size_t i = 0; i < filenames[0].size(); ++i ) {
				std::vector< wildcard::FilenameMatch > slice( filenames.size() ) ;
				for( std::size_t j = 0; j < slice.size(); ++j ) {
					slice[j] = filenames[j][i] ;
				}
				chain->add_source( std::auto_ptr< SNPDataSource >( SNPDataSourceRack::create( slice ))) ;
				if( notify_progress ) {
					notify_progress( ++count, total_number_of_files ) ;
				}
			}
		}
		return chain ;
	}
		
	SNPDataSourceChain::SNPDataSourceChain()
		: m_current_source(0), m_number_of_samples(0), m_moved_to_next_source_callback(0)
	{}

	SNPDataSourceChain::~SNPDataSourceChain() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	SNPDataSource::Metadata SNPDataSourceChain::get_metadata() const {
		return m_metadata ;
	}

	void SNPDataSourceChain::add_source( std::auto_ptr< SNPDataSource > source ) {
		if( m_sources.empty() ) {
			m_number_of_samples = source->number_of_samples() ;
			m_current_source = 0 ;
		}
		else if( source->number_of_samples() != m_number_of_samples ) {
			throw FileContainsSNPsOfDifferentSizes() ;
		}
		m_sources.push_back( source.release() ) ;
		//m_sources.back()->reset_to_start() ;

		// update metadata
		Metadata new_source_metadata = m_sources.back()->get_metadata() ;
		Metadata new_metadata ;
		std::set_union(
			m_metadata.begin(), m_metadata.end(),
			new_source_metadata.begin(), new_source_metadata.end(),
			std::inserter( new_metadata, new_metadata.end() )
		) ;
		m_metadata = new_metadata ;
	}

	unsigned int SNPDataSourceChain::number_of_samples() const {
		return m_number_of_samples ;
	}

	bool SNPDataSourceChain::has_sample_ids() const {
		return m_sources.size() > 0 && m_sources[0]->has_sample_ids() ;
	}

	void SNPDataSourceChain::get_sample_ids( GetSampleIds setter ) const {
		if( m_sources.size() > 0 ) {
			m_sources[0]->get_sample_ids( setter ) ;
		}
	}

	unsigned int SNPDataSourceChain::number_of_sources() const {
		return m_sources.size() ;
	}

	SNPDataSource::OptionalSnpCount SNPDataSourceChain::total_number_of_snps() const {
		unsigned int total_number_of_snps = 0 ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			if( !m_sources[i]->total_number_of_snps() ) {
				return OptionalSnpCount() ;
			}
			total_number_of_snps += *( m_sources[i]->total_number_of_snps() ) ;
		}
		return total_number_of_snps ;
	}

	SNPDataSource::OptionalSnpCount SNPDataSourceChain::number_of_snps_in_source( std::size_t source_index ) const {
		assert( source_index < m_sources.size() ) ;
		return m_sources[ source_index ]->total_number_of_snps() ;
	}

	SNPDataSource const& SNPDataSourceChain::get_source( std::size_t source_index ) const {
		assert( source_index < m_sources.size() ) ;
		return *m_sources[ source_index ] ;
	}

	SNPDataSourceChain::operator bool() const {
		if( m_current_source < m_sources.size() ) {	
			return *m_sources[ m_current_source ] ;
		}
		else {
			return false ;
		}
	}

	std::string SNPDataSourceChain::get_source_spec() const {
		std::string result = "chain:" ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			if( i > 0 ) {
				result += "," ;
			}
			result += m_sources[i]->get_source_spec() ;
		}
		return result ;
	}
	
	std::string SNPDataSourceChain::get_summary( std::string const& prefix, std::size_t width ) const {
		std::ostringstream ostr ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			ostr
				<< prefix << std::setw( width ) << "" << " (" ;
			if( number_of_snps_in_source( i ) ) {
				ostr << std::setw(7) << *number_of_snps_in_source( i ) << " snps" ;
			}
			else {
				ostr << "not computed" ;
			}
			ostr
				<< ")  "
				<< "\"" << m_sources[i]->get_source_spec()
				<< "\"\n" ;
		}
		if( total_number_of_snps() ) {
			ostr << prefix << std::setw( width ) << "" << " (total " << *total_number_of_snps() << " snps in " << m_sources.size() << " sources).\n" ;
		} else {
			ostr << prefix << std::setw( width ) << "" << " (total " << m_sources.size() << " sources, number of snps not computed).\n" ;
		}
		ostr << prefix << std::setw( width ) << "Number of samples:" << " " << number_of_samples() << "\n" ;
		return ostr.str() ;
	}
	

	void SNPDataSourceChain::reset_to_start_impl() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[i]->reset_to_start() ;
		}
		m_current_source = 0 ;
	}

	void SNPDataSourceChain::get_snp_identifying_data_impl( 
		VariantIdentifyingData* result
	) {
		move_to_next_nonempty_source_if_necessary() ;
		if( m_current_source < m_sources.size() ) {
			m_sources[m_current_source]->get_snp_identifying_data( result ) ;
		}
	}

	VariantDataReader::UniquePtr SNPDataSourceChain::read_variant_data_impl() {
		assert( m_current_source < m_sources.size() ) ;
		return m_sources[m_current_source]->read_variant_data() ;
	}

	void SNPDataSourceChain::ignore_snp_probability_data_impl(
	) {
		assert( m_current_source < m_sources.size() ) ;
		m_sources[m_current_source]->ignore_snp_probability_data() ;
	}

	void SNPDataSourceChain::set_moved_to_next_source_callback( moved_to_next_source_callback_t callback ) { m_moved_to_next_source_callback = callback ; }

	void SNPDataSourceChain::move_to_next_source() {
		++m_current_source ;
		if( m_moved_to_next_source_callback ) {
			m_moved_to_next_source_callback( m_current_source ) ;
		}
	}

	void SNPDataSourceChain::move_to_next_nonempty_source_if_necessary() {
		VariantIdentifyingData variant ;
		while(
			m_current_source < m_sources.size()
			&&
			!m_sources[ m_current_source ]->get_snp_identifying_data( &variant )
		) {
			move_to_next_source() ;
		}
	}
}

