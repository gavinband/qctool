
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <set>
#include <boost/unordered_map.hpp>
#include <boost/format.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/HLAIMPAsBiallelicVariantDataSource.hpp"
#include "genfile/string_utils/slice.hpp"

namespace genfile {
	HLAIMPAsBiallelicVariantDataSource::HLAIMPAsBiallelicVariantDataSource( std::string const& filename ):
		m_filename( filename ),
		m_allele_index( 0 )
	{
		setup( m_filename ) ;
	}

	HLAIMPAsBiallelicVariantDataSource::HLAIMPAsBiallelicVariantDataSource( std::istream& stream ):
		m_filename( "(unknown stream)" ),
		m_allele_index( 0 )
	{
		setup( stream ) ;
	}

	void HLAIMPAsBiallelicVariantDataSource::setup( std::string const& filename ) {
		std::auto_ptr< std::istream > stream_ptr = open_text_file_for_input( filename ) ;
		setup( *stream_ptr ) ;
	}

	void HLAIMPAsBiallelicVariantDataSource::setup( std::istream& stream ) {
		using string_utils::slice ;
		// Read column header
		std::string line ;
		std::getline( stream, line ) ;
		if( !stream ) {
			throw genfile::MalformedInputError(
				m_filename,
				"No lines in file",
				0
			) ;
		}

		std::vector< slice > const columns = slice( line ).strip( "\n\r" ).split( " \t" ) ;
		std::vector< std::string > alleles_in_order ;
		{
			std::vector< slice > const& elts = columns ;
			if( elts[0] != "IndividualID" ) {
				std::cerr << "!! " << elts[0] << "\n" ;
				throw genfile::MalformedInputError(
					m_filename,
					"Expected first column to be called \"IndividualID\"",
					0, 0
				) ;
			}
			// We now expect K*(K+1)/2 entries for all genotype combinations of K alleles
			// They are named allele1_allele2 so we can just split by underscore to get the allele names
			// We keep alleles in the order they appear in the input file.
			{
				std::map< std::string, std::size_t > alleles ;
				for( std::size_t i = 1; i < elts.size(); ++i ) {
					std::vector< slice > theseAlleles = elts[i].split( "_" ) ;
					std::string const allele1 = theseAlleles[0] ;
					if( alleles.find( allele1 ) == alleles.end() ) {
						alleles.insert( std::make_pair( allele1, alleles_in_order.size() ) ) ;
						alleles_in_order.push_back( allele1 ) ;
					}
				}
				assert( alleles.size() == alleles_in_order.size() ) ;
			}
			std::size_t const K = alleles_in_order.size() ;
			
			// The following check is not used as files don't actually contain all allele combinations.
			std::size_t const expectedColumnCount = 1 + ( K*(K+1)/2 ) ;
			if( elts.size() != expectedColumnCount ) {
				std::cerr << "!! " << ( boost::format( "Unexpected number of columns (%d, expected %d*(%d+1)/2+1 = %d)" ) % elts.size() % K % K % expectedColumnCount ) << "\n" ;
//				throw genfile::MalformedInputError(
//					m_filename,
//					( boost::format( "Wrong number of columns (%d, expected %d*(%d+1)/2+1 = %d)" ) % elts.size() % K % K % expectedColumnCount ).str(),
//					0
//				) ;
			}
			m_alleles.swap( alleles_in_order ) ;
			
			// Now compute the indices of homs and hets for each allele
			m_allele_dosage_columns.resize( m_alleles.size() ) ;
			for( std::size_t allele_i = 0; allele_i < m_alleles.size(); ++allele_i ) {
				std::string const& allele = m_alleles[ allele_i ] ;
				m_allele_dosage_columns[ allele_i ].resize( 3 ) ;
				for( std::size_t i = 1; i < columns.size(); ++i ) {
					if( columns[i] == allele + "_" + allele ) {
						m_allele_dosage_columns[ allele_i ][2].push_back( i-1 ) ;
					} else if(
						(columns[i].substr( 0, allele.size() + 1 ) == ( allele + "_" ))
						|| ( columns[i].size() > ( allele.size() + 1 ) && columns[i].substr( columns[i].size() - allele.size() - 1, columns[i].size() ) == ( "_" + allele ))
					) {
						m_allele_dosage_columns[ allele_i ][1].push_back( i-1 ) ;
					} else {
						m_allele_dosage_columns[ allele_i ][0].push_back( i-1 ) ;
					}
				}
			}
		}
		
		// Read sample IDs and data.
		{
			double const NA = std::numeric_limits< double >::quiet_NaN() ;
			std::vector< Eigen::VectorXd > data ;
			std::vector< std::string > samples ;
			std::set< std::string > samplesSet ;
			
			for( std::size_t lineCount = 1; std::getline( stream, line ); ++lineCount ) {
				std::vector< slice > const elts = slice( line ).strip( "\n\r" ).split( " \t" ) ;
				if( elts.size() != columns.size() ) {
					throw genfile::MalformedInputError(
						m_filename,
						( boost::format( "Wrong number of entries (%d, expected %d)" ) % elts.size() % columns.size() ).str(),
						lineCount
					) ;
				}
			
				if( samplesSet.find( elts[0] ) != samplesSet.end() ) {
					throw genfile::MalformedInputError(
						m_filename,
						( boost::format( "Expected unique sample IDs (sample \"%s\" is repeated)" ) % std::string( elts[0] ) ).str(),
						lineCount, 0
					) ;
				}
				samples.push_back( elts[0] ) ;
				samplesSet.insert( elts[0] ) ;

				data.push_back( Eigen::VectorXd::Constant( elts.size() - 1 , NA )) ;
				for( std::size_t i = 1; i < elts.size(); ++i ) {
					if( elts[i] != "NA" ) {
						try {
							data.back()(i-1) = genfile::string_utils::to_repr< double >( elts[i] ) ;
						} catch( genfile::string_utils::StringConversionError const& ) {
							throw genfile::MalformedInputError(
								m_filename,
								( boost::format( "Expected a numerical value, not \"%s\"" ) % std::string( elts[i] ) ).str(),
								lineCount, i
							) ;
						}
					}
				}
			}
			m_samples.swap( samples ) ;
			m_data.swap( data ) ;
		}
	}

	SNPDataSource::Metadata HLAIMPAsBiallelicVariantDataSource::get_metadata() const {
		SNPDataSource::Metadata result ;

		{
			std::map< std::string, std::string > format ;
			format[ "ID" ] = "GP" ;
			format[ "Type" ] = "Float" ;
			format[ "Number" ] = "G" ;
			format[ "Description" ] = "Genotype probability" ;
			result.insert( std::make_pair( "FORMAT", format )) ;
		}

		return result ;
	}

	unsigned int HLAIMPAsBiallelicVariantDataSource::number_of_samples() const {
		return m_samples.size() ;
	}

	void HLAIMPAsBiallelicVariantDataSource::get_sample_ids( GetSampleIds getter ) const {
		for( std::size_t i = 0; i < m_samples.size(); ++i ) {
			getter( i, m_samples[i] ) ;
		}
	}

	HLAIMPAsBiallelicVariantDataSource::OptionalSnpCount HLAIMPAsBiallelicVariantDataSource::total_number_of_snps() const {
		return m_alleles.size() ;
	}

	HLAIMPAsBiallelicVariantDataSource::operator bool() const {
		return m_allele_index < m_alleles.size() ;
	}

	std::string HLAIMPAsBiallelicVariantDataSource::get_source_spec() const { return "HLAIMPAsBiallelicVariantDataSource(\"" + m_filename + "\")" ; }

	void HLAIMPAsBiallelicVariantDataSource::reset_to_start_impl() {
		m_exhausted = false ;
	}
	
	void HLAIMPAsBiallelicVariantDataSource::read_snp_identifying_data_impl( 
		uint32_t* number_of_samples,
		std::string* SNPID,
		std::string* RSID,
		Chromosome* chromosome,
		uint32_t* SNP_position,
		std::string* allele1,
		std::string* allele2
	) {
		if( m_allele_index < m_alleles.size() ) {
			*number_of_samples = m_samples.size() ;
			*SNPID = m_alleles[ m_allele_index ] ;
			*RSID = m_alleles[ m_allele_index ] ;
			*chromosome = Chromosome( "06" ) ;
			*SNP_position = 0 ;
			*allele1 = "non-" + m_alleles[ m_allele_index ] + "_alleles" ;
			*allele2 = m_alleles[ m_allele_index ] ;
		}
	}

	namespace impl {
		struct HLAIMPVariantDataReader: public VariantDataReader {
			static UniquePtr create( HLAIMPAsBiallelicVariantDataSource const& source, std::size_t allele_index ) {
				return UniquePtr( new HLAIMPVariantDataReader( source, allele_index )) ;
			}
			
			HLAIMPVariantDataReader( HLAIMPAsBiallelicVariantDataSource const& source, std::size_t allele_index ):
				m_source( source ),
				m_allele_index( allele_index ),
				m_dosage_columns( m_source.m_allele_dosage_columns[ m_allele_index ] )
			{}
			
			HLAIMPVariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				std::size_t const K = m_source.m_alleles.size() ;
				std::size_t const genotypeCount = K*(K+1)/2 ;
				if( spec == ":genotypes:" || spec == "GP" ) {
					setter.initialise( m_source.m_samples.size(), 2 ) ;
					for( std::size_t sample_i = 0; sample_i < m_source.m_samples.size(); ++sample_i ) {
						if( setter.set_sample( sample_i )) {
							setter.set_number_of_entries( 3, ePerUnorderedGenotype, eProbability ) ;
							for( std::size_t g = 0; g < 3; ++g ) {
								double value = 0 ;
								for( std::size_t j = 0; j < m_dosage_columns[g].size(); ++j ) {
									value += m_source.m_data[sample_i](m_dosage_columns[g][j]) ;
								}
								setter.set_value( value ) ;
							}
						}
					}
				} else {
					throw BadArgumentError(
						"genfile::HLAIMPVariantDataReader::get()",
						"spec=\"" + spec + "\"",
						"Only \"GP\" and \":genotypes:\" are supported in HLAIMP data."
					) ;
				}
				return *this ;
			}

			std::size_t get_number_of_samples() const {
				return m_source.m_samples.size() ;
			}
		
			bool supports( std::string const& spec ) const {
				return spec == "GP" || spec == ":genotypes:" ;
			}
		
			void get_supported_specs( SpecSetter setter ) const {
				setter( "GP", "Float" ) ;
				setter( ":genotypes:", "Float" ) ;
			}
		private:
			HLAIMPAsBiallelicVariantDataSource const& m_source ;
			std::size_t const m_allele_index ;
			std::vector< std::vector< std::size_t > > const& m_dosage_columns ;
		} ;
	}

	VariantDataReader::UniquePtr HLAIMPAsBiallelicVariantDataSource::read_variant_data_impl() {
		assert( !m_exhausted ) ;
		return impl::HLAIMPVariantDataReader::create( *this, m_allele_index++ ) ;
	}

	void HLAIMPAsBiallelicVariantDataSource::ignore_snp_probability_data_impl() {
		++m_allele_index ;
	}
}

