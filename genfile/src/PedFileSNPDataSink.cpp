#include <vector>
#include <utility>
#include <fstream>
#include "genfile/SNPDataSink.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/Pedigree.hpp"
#include "genfile/PedFileSNPDataSink.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	namespace impl {
		char get_representation_of_allele( char allele ) {
			switch( allele ) {
				case 'A': return '1'; break ;
				case 'C': return '2'; break ;
				case 'G': return '3'; break ;
				case 'T': return '4'; break ;
			}
			return '?' ;
		}
		
		std::string sex_to_string( Pedigree::Sex sex ) {
			switch( sex ) {
				case Pedigree::eMale: return "M"; break ;
				case Pedigree::eFemale: return "F"; break ;
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

	PedFileSNPDataSink::PedFileSNPDataSink(
		CohortIndividualSource const& samples,
		Pedigree const& pedigree,
		std::string const& output_filename,
		double call_threshhold
	):
		m_samples( samples ),
		m_pedigree( pedigree ),
		m_phenotypes( impl::get_phenotypes( samples )),
		m_pedigree_to_sample_mapping( get_pedigree_to_sample_mapping( pedigree, samples )),
		m_output_filename( output_filename ),
		m_call_threshhold( call_threshhold )
	{
		assert( call_threshhold > 0.5 ) ;
		if( m_output_filename.size() >= 4 &&
			(
				( m_output_filename.substr( m_output_filename.size() - 4, 4 ) == ".ped" )
				||
				( m_output_filename.substr( m_output_filename.size() - 4, 4 ) == ".map" )
			)
		) {
			throw BadArgumentError( "PedFileSNPDataSink::PedFileSNPDataSink", "output_filename = \"" + output_filename + "\"" ) ;
		}
	}
	
	std::map< std::string, std::size_t > PedFileSNPDataSink::get_pedigree_to_sample_mapping(
		Pedigree const& pedigree,
		CohortIndividualSource const& samples
	) {
		std::map< std::string, std::size_t > result ;

		std::vector< std::string > ID_1( samples.get_number_of_individuals() ) ;
		std::vector< std::string > ID_2( samples.get_number_of_individuals() ) ;
		
		for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
			ID_1[i] = samples.get_entry( i, "ID_1" ).as< std::string >() ;
			ID_2[i] = samples.get_entry( i, "ID_2" ).as< std::string >() ;
		}
		
		for( std::size_t i = 0; i < pedigree.get_number_of_individuals(); ++i ) {
			std::string const& id = pedigree.get_id_of( i ) ;
			//std::cerr << "Matching up individual " << i << ": " << id << ".\n" ;
			std::vector< std::string > const* which_id = &ID_1 ;
			std::vector< std::string >::const_iterator where = std::find( which_id->begin(), which_id->end(), id ) ;
			if( where == ID_1.end() ) {
				// look in ID 2 instead
				which_id = &ID_2 ;
				where = std::find( which_id->begin(), which_id->end(), id ) ;
			}
			if( where != which_id->end() ) {
				// check match is unique.
				std::vector< std::string >::const_iterator next = where ;
				++next ;
				assert( std::find( next, which_id->end(), id ) == which_id->end() ) ;
				result[ id ] = ( where - which_id->begin() )  ;
			}
		}
		return result ;
	}
	
	PedFileSNPDataSink::~PedFileSNPDataSink() {
		write_ped_file( m_output_filename + ".ped" ) ;
		write_map_file( m_output_filename + ".map" ) ;
	}
	
	void PedFileSNPDataSink::write_ped_file( std::string const& output_filename ) const {
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
				<< " " << impl::sex_to_string( m_pedigree.get_sex_of( id ) ) ;

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
	
	void PedFileSNPDataSink::write_map_file( std::string const& output_filename ) const {
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
	
	void PedFileSNPDataSink::write_snp_impl(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		Chromosome chromosome,
		uint32_t SNP_position,
		char first_allele,
		char second_allele,
		GenotypeProbabilityGetter const& get_AA_probability,
		GenotypeProbabilityGetter const& get_AB_probability,
		GenotypeProbabilityGetter const& get_BB_probability
	) {
		assert( m_written_snps.size() == m_written_alleles.size() ) ;
		assert( number_of_samples == m_samples.get_number_of_individuals() ) ;
		
		m_written_snps.push_back(
			SNPIdentifyingData(
				SNPID,
				RSID,
				GenomePosition( chromosome, SNP_position ),
				first_allele,
				second_allele
			)
		) ;

		SNPIdentifyingData const& snp = m_written_snps.back() ;

		m_written_alleles.push_back( std::vector< std::pair< char, char > >() ) ;
		// Convert probabilities to calls.
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			double AA = get_AA_probability( i ),
				AB = get_AB_probability( i ),
				BB = get_BB_probability( i ) ;
			std::pair< char, char > alleles ;
			// use the folloiwn
			if( AA > m_call_threshhold ) {
				alleles.first = alleles.second = impl::get_representation_of_allele( snp.get_first_allele() ) ;
			}
			else if( AB > m_call_threshhold ) {
				alleles.first = impl::get_representation_of_allele( snp.get_first_allele() ) ;
				alleles.second = impl::get_representation_of_allele( snp.get_second_allele() ) ;
			}
			else if( BB > m_call_threshhold ) {
				alleles.first = alleles.second = impl::get_representation_of_allele( snp.get_second_allele() ) ;
			}
			else {
				alleles.first = alleles.second = '?' ;
			}
			m_written_alleles.back().push_back( alleles ) ;
		}
		assert( m_written_snps.size() == m_written_alleles.size() ) ;
		assert( m_written_alleles.back().size() == m_samples.get_number_of_individuals() ) ;
	}
}