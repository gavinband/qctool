
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <utility>
#include <fstream>
#include "genfile/SNPDataSink.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/Pedigree.hpp"
#include "genfile/PedigreeMappingPedFileSNPDataSink.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	namespace {
		std::string sex_to_string( Pedigree::Sex sex ) {
			switch( sex ) {
				case Pedigree::eMale: return "1"; break ;
				case Pedigree::eFemale: return "2"; break ;
				case Pedigree::eUnknown: return "NA"; break ;
			}
			return "NA";
		}
		
		std::vector< std::string > get_phenotypes( genfile::CohortIndividualSource const& samples ) {
			std::vector< std::string > result ;
			genfile::CohortIndividualSource::ColumnSpec column_spec = samples.get_column_spec() ;
			for( std::size_t i = 0; i < column_spec.size(); ++i ) {
				if( column_spec.get_spec(i).is_phenotype() ) {
					result.push_back( column_spec[i].name() ) ;
				}
			}
			return result ;
		}
	}

	PedigreeMappingPedFileSNPDataSink::PedigreeMappingPedFileSNPDataSink(
		CohortIndividualSource const& samples,
		Pedigree const& pedigree,
		std::string const& output_ped_filename,
		double call_threshhold
	):
		m_samples( samples ),
		m_pedigree( pedigree ),
		m_phenotypes( get_phenotypes( samples )),
		m_pedigree_to_sample_mapping( get_pedigree_to_sample_mapping( pedigree, samples )),
		m_call_threshhold( call_threshhold )
	{
		assert( call_threshhold > 0.5 ) ;
		if( output_ped_filename.size() < 4 || output_ped_filename.substr( output_ped_filename.size() - 4, 4 ) != ".ped" ) {
			throw BadArgumentError( "PedigreeMappingPedFileSNPDataSink::PedigreeMappingPedFileSNPDataSink", "output_ped_filename = \"" + output_ped_filename + "\"" ) ;
		}
		if( m_pedigree_to_sample_mapping.empty() ) {
			// No matches between sample file and pedigree could be found.
			// All output genotypes would be NA, let's throw an error in this case.
			throw MismatchError( "genfile::PedigreeMappingPedFileSNPDataSink::PedigreeMappingPedFileSNPDataSink()", "Pedigree / sample file", "Individual id (column 2 of pedigree file)", "ID_1 / ID_2 column." ) ;
		}
		m_output_filename_stub = output_ped_filename.substr( 0, output_ped_filename.size() - 4 ) ;
	}
	
	std::string PedigreeMappingPedFileSNPDataSink::get_spec() const {
		return m_output_filename_stub + "[.ped|.dat|.map]";
	}
	
	std::map< std::string, std::size_t > PedigreeMappingPedFileSNPDataSink::get_pedigree_to_sample_mapping(
		Pedigree const& pedigree,
		CohortIndividualSource const& samples
	) {
		std::map< std::string, std::size_t > result ;

		std::vector< std::string > ID_1( samples.get_number_of_individuals() ) ;
		
		for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
			ID_1[i] = samples.get_entry( i, "ID_1" ).as< std::string >() ;
		}
		
		for( std::size_t i = 0; i < pedigree.get_number_of_individuals(); ++i ) {
			std::string const& id = pedigree.get_id_of( i ) ;
			//std::cerr << "Matching up individual " << i << ": " << id << ".\n" ;
			std::vector< std::string >::const_iterator where = std::find( ID_1.begin(), ID_1.end(), id ) ;
			if( where !=ID_1.end() ) {
				// check match is unique.
				std::vector< std::string >::const_iterator next = where ;
				std::vector< std::string >::const_iterator const end = ID_1.end() ;
				++next ;
				assert( std::find( next, end, id ) == end ) ;
				result[ id ] = ( where - ID_1.begin() )  ;
			}
		}
		return result ;
	}
	
	PedigreeMappingPedFileSNPDataSink::~PedigreeMappingPedFileSNPDataSink() {
		write_ped_file( m_output_filename_stub + ".ped" ) ;
		write_dat_file( m_output_filename_stub + ".dat" ) ;
		write_map_file( m_output_filename_stub + ".map" ) ;
	}
	
	void PedigreeMappingPedFileSNPDataSink::write_ped_file( std::string const& output_filename ) const {
		assert( m_written_snps.size() == m_written_alleles.size() ) ;
		for( std::size_t i = 0; i < m_written_alleles.size(); ++i ) {
			assert( m_written_alleles[i].size() == m_samples.get_number_of_individuals() ) ;
		}

		std::auto_ptr< std::ostream> out(
			open_text_file_for_output(
				output_filename,
				get_compression_type_indicated_by_filename( output_filename )
			)
		) ;
		for( std::size_t i = 0; i < m_pedigree.get_number_of_individuals(); ++i ) {
			// first five columns are:
			// family, id, father, mother, sex
			std::string const id = m_pedigree.get_id_of( i ) ;
			(*out) << m_pedigree.get_family_of( id )
				<< " " << id
				<< " " << m_pedigree.get_parents_of( id )[0]
				<< " " << m_pedigree.get_parents_of( id )[1]
				<< " " << sex_to_string( m_pedigree.get_sex_of( id ) ) ;

			std::map< std::string, std::size_t >::const_iterator sample = m_pedigree_to_sample_mapping.find( id ) ;
			if( sample != m_pedigree_to_sample_mapping.end() ) {
				// Next columns are phenotypes
				for( std::size_t phenotype_i = 0; phenotype_i < m_phenotypes.size(); ++phenotype_i ) {
					(*out) << " " << m_samples.get_entry( sample->second, m_phenotypes[ phenotype_i ] ) ;
				}
				// Next columns are the SNPs, which we go through in order.
				for( std::size_t i = 0; i < m_written_alleles.size(); ++i ) {
					(*out) << " "
						<< m_written_alleles[ i ][ sample->second ].first
						<< "/"
						<< m_written_alleles[ i ][ sample->second ].second ;
				}
			}
			else {
				for( std::size_t phenotype_i = 0; phenotype_i < m_phenotypes.size(); ++phenotype_i ) {
					(*out) << " NA" ;
				}
				for( std::size_t i = 0 ; i < m_written_alleles.size(); ++i ) {
					(*out) << " NA/NA" ;
				}
			}
			(*out) << "\n" ;
		}
		assert( (*out) ) ;
	}
	
	void PedigreeMappingPedFileSNPDataSink::write_dat_file( std::string const& output_filename ) const {
		std::auto_ptr< std::ostream> out(
			open_text_file_for_output(
				output_filename,
				get_compression_type_indicated_by_filename( output_filename )
			)
		) ;
		for( std::size_t phenotype_i = 0; phenotype_i < m_phenotypes.size(); ++phenotype_i ) {
			(*out) << "T " << m_phenotypes[ phenotype_i ] << "\n" ;
		}
		for( std::size_t i = 0 ; i < m_written_snps.size(); ++i ) {
			(*out) << "M " << m_written_snps[i].get_rsid() << "\n" ;
		}
		assert( (*out) ) ;
	}

	void PedigreeMappingPedFileSNPDataSink::write_map_file( std::string const& output_filename ) const {
		std::auto_ptr< std::ostream > file = open_text_file_for_output( output_filename ) ;
		for( std::size_t i = 0; i < m_written_snps.size(); ++i ) {
			VariantIdentifyingData const& snp = m_written_snps[i] ;
			(*file)
				<< snp.get_position().chromosome()
				<< " "
				<< snp.get_rsid()
				<< " "
				<< snp.get_position().position()
				<< "\n" ;
		}
	}
	
	void PedigreeMappingPedFileSNPDataSink::write_snp_impl(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		Chromosome chromosome,
		uint32_t SNP_position,
		std::string first_allele,
		std::string second_allele,
		GenotypeProbabilityGetter const& get_AA_probability,
		GenotypeProbabilityGetter const& get_AB_probability,
		GenotypeProbabilityGetter const& get_BB_probability,
		Info const& info
	) {
		assert( m_written_snps.size() == m_written_alleles.size() ) ;
		assert( number_of_samples == m_samples.get_number_of_individuals() ) ;
		
		m_written_snps.push_back(
			VariantIdentifyingData(
				SNPID,
				RSID,
				GenomePosition( chromosome, SNP_position ),
				first_allele,
				second_allele
			)
		) ;

		m_written_alleles.push_back( std::vector< std::pair< char, char > >() ) ;
		// Convert probabilities to calls.
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			double AA = get_AA_probability( i ),
				AB = get_AB_probability( i ),
				BB = get_BB_probability( i ) ;
			std::pair< char, char > alleles = std::make_pair( '?', '?' );
			if( AA > m_call_threshhold ) {
				alleles.first = alleles.second = '1' ;
			}
			else if( AB > m_call_threshhold ) {
				alleles.first = '1' ;
				alleles.second = '2' ;
			}
			else if( BB > m_call_threshhold ) {
				alleles.first = alleles.second = '2' ;
			}
			else {
				alleles.first = alleles.second = 'x' ;
			}
			m_written_alleles.back().push_back( alleles ) ;
		}
		assert( m_written_snps.size() == m_written_alleles.size() ) ;
		assert( m_written_alleles.back().size() == m_samples.get_number_of_individuals() ) ;
	}
}
