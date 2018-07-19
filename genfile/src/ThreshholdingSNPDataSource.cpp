
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <queue>
#include <memory>
#include <boost/format.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/ThreshholdingSNPDataSource.hpp"
#include "genfile/GPThresholdingGTSetter.hpp"

namespace genfile {
// This SNPDataSource applies a threshhold to genotype calls.
	ThreshholdingSNPDataSource::ThreshholdingSNPDataSource( SNPDataSource::UniquePtr source, double const threshhold ):
		m_source( source ),
		m_threshhold( threshhold )
	{}
	
	ThreshholdingSNPDataSource::operator bool() const {
		return *m_source ;
	}

	SNPDataSource::Metadata ThreshholdingSNPDataSource::get_metadata() const {
		// Remove GT from metadata
		Metadata metadata = m_source->get_metadata() ;
		for( Metadata::iterator i = metadata.begin(); i != metadata.end(); ) {
			Metadata::iterator erase_i = i++ ;
			if( erase_i->first == "FORMAT" ) {
				if( erase_i->second.find( "ID" ) != erase_i->second.end() && erase_i->second[ "ID" ] == "GT" ) {
					metadata.erase( erase_i ) ;
				}
			}
		}
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "GT" ;
		format[ "Number" ] = "1" ;
		format[ "Type" ] = "String" ;
		format[ "Description" ] = ( boost::format( "Genotype call probabilities, threshholded at %.2f" ) % m_threshhold ).str() ;
		metadata.insert( std::make_pair( "FORMAT", format )) ;
		return metadata ;
	}

	unsigned int ThreshholdingSNPDataSource::number_of_samples() const {
		return m_source->number_of_samples() ;
	}
	void ThreshholdingSNPDataSource::get_sample_ids( GetSampleIds getter ) const {
		return m_source->get_sample_ids( getter ) ;
	}
	SNPDataSource::OptionalSnpCount ThreshholdingSNPDataSource::total_number_of_snps() const {
		return m_source->total_number_of_snps() ;
	}
	std::string ThreshholdingSNPDataSource::get_source_spec() const {
		return "ThreshholdingSNPDataSource(" + m_source->get_source_spec() + ")" ;
	}
	SNPDataSource const& ThreshholdingSNPDataSource::get_parent_source() const {
		return *m_source ;
	}
	SNPDataSource const& ThreshholdingSNPDataSource::get_base_source() const {
		return m_source->get_base_source() ;
	}
	std::string ThreshholdingSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return "ThreshholdingSNPDataSource(" + m_source->get_source_spec() + ")" ;
	}
	void ThreshholdingSNPDataSource::get_snp_identifying_data_impl( VariantIdentifyingData* variant	) {
		m_source->get_snp_identifying_data( variant ) ;
	}

	namespace {
		struct ThreshholdingDataReader: public VariantDataReader {
			ThreshholdingDataReader( VariantDataReader::UniquePtr source, double const threshhold ):
				m_source( source ),
				m_threshhold( threshhold )
			{}

			~ThreshholdingDataReader() {}

			typedef GPThresholdingGTSetter Threshholder ;

			VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				if( spec == ":genotypes:" || spec == "GT" ) {
					if( !m_source->supports( "GP" )) {
						throw genfile::BadArgumentError(
							"genfile::ThreshholdingDataReader::Threshholder::get()",
							"setter",
							"Underlying source must support genotype probabilities (GP) field."
						) ;
					}
					Threshholder threshholder( setter, m_threshhold ) ;
					m_source->get( "GP", threshholder ) ;
				} else {
					m_source->get( spec, setter ) ;
				}
				return *this ;
			}

			bool supports( std::string const& spec ) const {
				return spec == ":genotypes:" || spec == "GT" || m_source->supports( spec ) ;
			}
			void get_supported_specs( SpecSetter setter ) const {
				setter( ":genotypes:", "Integer" ) ;
				setter( "GT", "Integer" ) ;
				return m_source->get_supported_specs( setter ) ;
			}
			std::size_t get_number_of_samples() const {
				return m_source->get_number_of_samples() ;
			}
		private:
			VariantDataReader::UniquePtr m_source ;
			std::string const m_field ;
			double const m_threshhold ;
		} ;
	}

	VariantDataReader::UniquePtr ThreshholdingSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr(
			new ThreshholdingDataReader(
				m_source->read_variant_data(),
				m_threshhold
			)
		) ;
	}

	void ThreshholdingSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	void ThreshholdingSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}

