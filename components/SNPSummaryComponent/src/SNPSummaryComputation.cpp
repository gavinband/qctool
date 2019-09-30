
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "genfile/Error.hpp"
#include "genfile/ToGP.hpp"
#include "genfile/VariantDataReader.hpp"
#include "metro/likelihood/Multinomial.hpp"
#include "components/SNPSummaryComponent/SNPHWE.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/HWEComputation.hpp"
#include "components/SNPSummaryComponent/IntensitySummaryComputation.hpp"
#include "components/SNPSummaryComponent/ClusterFitComputation.hpp"
#include "components/SNPSummaryComponent/InfoComputation.hpp"

// #define DEBUG_SNP_SUMMARY_COMPUTATION 1

namespace stats {	
	namespace {
		struct AlleleCountClient: public genfile::VariantDataReader::PerSampleSetter {
			AlleleCountClient( std::vector< double >* counts ):
				m_counts( counts ),
				m_number_of_alleles( 0 )
			{
				assert( counts != 0 ) ;
			}

			~AlleleCountClient() throw() {}
			
			void set_counts( std::vector< double >* counts ) {
				assert( counts != 0 ) ;
				m_counts = counts ;
			}
			
			void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
				m_number_of_alleles = number_of_alleles ;
				assert( m_counts->size() == number_of_alleles ) ;
			}
			
			bool set_sample( std::size_t i ) {
				m_ploidy = 0 ;
				m_table = 0 ;
				return true ;
			}
			
			void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				assert( order_type == genfile::ePerUnorderedGenotype ) ;
				assert( value_type == genfile::eProbability ) ;
				
				std::map< uint32_t, genfile::impl::Enumeration >::iterator where = m_tables.find( ploidy ) ;
				if( where == m_tables.end() ) {
					std::pair< std::map< uint32_t, genfile::impl::Enumeration >::iterator, bool >
						result = m_tables.insert( std::make_pair( ploidy, genfile::impl::enumerate_unphased_genotypes( ploidy ) )) ;
					assert( result.second ) ;
					where = result.first ;
				}
				m_ploidy = ploidy ;
				m_table = &(where->second) ;
			}
			
			void set_value( std::size_t value_i, genfile::MissingValue const value ) {
				set_value( value_i, 0.0 ) ;
			}

			void set_value( std::size_t value_i, double const value ) {
				assert( value_i <= std::size_t( std::numeric_limits< uint16_t >::max() )) ;
				genfile::impl::Enumeration const& enumeration = *m_table ;
				std::size_t const& maxAlleles = enumeration.first.second ;
				assert( m_number_of_alleles <= maxAlleles ) ;
				uint16_t const encodedGenotype = enumeration.second.second[ uint16_t( value_i ) ] ;
				uint32_t const bitsPerAllele = enumeration.first.first ;
				uint16_t mask = uint16_t( 0xFFFF ) >> ( 16 - bitsPerAllele ) ;
				// Dosage of the 1st allele in the genotype is not encoded directly.
				// We compute it as the ploidy minus the dosage of other alleles.
				uint16_t dosage_of_nonref_alleles = 0 ;
				for( std::size_t allele = 1; allele < m_number_of_alleles; ++allele ) {
					uint16_t const dosage = ( encodedGenotype >> ( (allele-1) * bitsPerAllele )) & mask ;
					(*m_counts)[allele] += value * dosage ;
					dosage_of_nonref_alleles += dosage ;
				}
				(*m_counts)[0] += value * (m_ploidy - dosage_of_nonref_alleles) ;
			}
			
			void finalise() {}
		private:
			std::vector< double >* m_counts ;
			std::size_t m_number_of_alleles ;
			uint32_t m_ploidy ;
			std::map< uint32_t, genfile::impl::Enumeration > m_tables ;
			genfile::impl::Enumeration* m_table ;
		} ;
	}

	struct AlleleCountComputation: public SNPSummaryComputation {
	public:
		
		AlleleCountComputation():
			m_counter( &m_counts )
		{}

		void list_variables( NameCallback callback ) const {
			using genfile::string_utils::to_string ;
			// By default we list 10 alleles
			callback( "number_of_alleles" ) ;
			for( std::size_t i = 0; i < 10; ++i ) {
				callback( "allele" + to_string(i+1) + "_count" ) ;
			}
		}

		void operator()(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy,
			genfile::VariantDataReader& data_reader,
			ResultCallback callback
		) {
			using genfile::string_utils::to_string ;
			m_counts = std::vector< double >( snp.number_of_alleles(), 0.0 ) ;
			m_counter.set_counts( &m_counts ) ;
			callback( "number_of_alleles", int64_t( snp.number_of_alleles() )) ;
			data_reader.get( ":genotypes:", genfile::to_GP_unphased( m_counter ) ) ;
			std::size_t i = 0 ;
			for( ; i < snp.number_of_alleles(); ++i ) {
				callback( "allele" + to_string(i+1) + "_count", m_counts[i] ) ;
			}
			for( ; i < 10; ++i ) {
				callback( "allele" + to_string(i+1) + "_count", genfile::MissingValue() ) ;
			}
		}
		
		std::string get_summary( std::string const& prefix, std::size_t column_width ) const {
			return prefix + "AlleleCountComputation" ;
		}
	private:
		std::vector< double > m_counts ;
		AlleleCountClient m_counter ;
	} ;

	struct AlleleFrequencyComputation: public SNPSummaryComputation
	{
		AlleleFrequencyComputation( std::string const& what ):
			m_compute_counts( what == "everything" || what == "counts" ),
			m_compute_frequencies( what == "everything" )
		{
			assert( what == "counts" || what == "everything" ) ;
		}
		void operator()(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy,
			genfile::VariantDataReader&,
			ResultCallback callback
		) {
			if( snp.number_of_alleles() == 2 ) {
				bool const allDiploid = ( ploidy.array() == 2 ).cast< int >().sum() == ploidy.size() ;
				if( allDiploid ) {
					compute_autosomal_frequency( snp, genotypes, callback ) ;
				} else {
					compute_sex_chromosome_frequency( snp, genotypes, ploidy, callback ) ;
				}
			}
		}
		
		void compute_sex_chromosome_frequency(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy,
			ResultCallback callback
		) {
			Genotypes haploid_genotypes = genotypes ;
			Genotypes diploid_genotypes = genotypes ;
			assert( std::size_t( genotypes.rows() ) == ploidy.size() ) ;

			for( int i = 0; i < genotypes.rows(); ++i ) {
				if( ploidy(i) != 1 ) {
					haploid_genotypes.row(i).setZero() ;
				}
				if( ploidy(i) != 2 ) {
					diploid_genotypes.row(i).setZero() ;
				}
			}
			
			double const a_allele_count = haploid_genotypes.col(0).sum()
				+ ( ( 2.0 * diploid_genotypes.col(0).sum() ) + diploid_genotypes.col(1).sum() ) ;
			double const b_allele_count = haploid_genotypes.col(1).sum()
				+ ( ( 2.0 * diploid_genotypes.col(2).sum() ) + diploid_genotypes.col(1).sum() ) ;

			if( m_compute_counts ) {
				callback( "alleleA_count", a_allele_count ) ;
				callback( "alleleB_count", b_allele_count ) ;
			}

			if( m_compute_frequencies ) {
				double const total_allele_count = ( haploid_genotypes.sum() + 2.0 * diploid_genotypes.sum() ) ;
				double const a_allele_freq = a_allele_count / total_allele_count ;
				double const b_allele_freq = b_allele_count / total_allele_count ;

				callback( "alleleA_frequency", a_allele_freq ) ;
				callback( "alleleB_frequency", b_allele_freq ) ;

				if( a_allele_freq < b_allele_freq ) {
					callback( "minor_allele_frequency", a_allele_freq ) ;
					callback( "minor_allele", snp.get_allele(0) ) ;
					callback( "major_allele", snp.get_allele(1) ) ;
				}
				else if( a_allele_freq > b_allele_freq ) {
					callback( "minor_allele_frequency", b_allele_freq ) ;
					callback( "minor_allele", snp.get_allele(1) ) ;
					callback( "major_allele", snp.get_allele(0) ) ;
				} else {
					callback( "minor_allele_frequency", a_allele_freq ) ;
				}
			}
		}

		void compute_autosomal_frequency( VariantIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) {
			double const a_allele_count = ( 2.0 * genotypes.col(0).sum() ) + genotypes.col(1).sum() ;
			double const b_allele_count = ( 2.0 * genotypes.col(2).sum() ) + genotypes.col(1).sum() ;

			if( m_compute_counts ) {
				callback( "alleleA_count", a_allele_count ) ;
				callback( "alleleB_count", b_allele_count ) ;
			}
			
			if( m_compute_frequencies ) {
				double const total_allele_count = ( 2.0 * genotypes.sum() ) ;
				double const a_allele_freq = a_allele_count / total_allele_count ;
				double const b_allele_freq = b_allele_count / total_allele_count ;

				callback( "alleleA_frequency", a_allele_freq ) ;
				callback( "alleleB_frequency", b_allele_freq ) ;

				if( a_allele_freq < b_allele_freq ) {
					callback( "minor_allele_frequency", a_allele_freq ) ;
					callback( "minor_allele", snp.get_allele(0) ) ;
					callback( "major_allele", snp.get_allele(1) ) ;
				}
				else if( a_allele_freq > b_allele_freq ) {
					callback( "minor_allele_frequency", b_allele_freq ) ;
					callback( "minor_allele", snp.get_allele(1) ) ;
					callback( "major_allele", snp.get_allele(0) ) ;
				} else {
					callback( "minor_allele_frequency", a_allele_freq ) ;
				}
			}
		}
		
		std::string get_summary( std::string const& prefix, std::size_t column_width ) const {
			return prefix + "AlleleFrequencyComputation" ;
		}
	private:
		bool const m_compute_counts ;
		bool const m_compute_frequencies ;
	} ;

	// What proportion of the mass on a genotype is due to high-confidence calls?
	struct CallMassComputation: public SNPSummaryComputation
	{
		
		CallMassComputation( double const threshhold = 0.9 ):
			m_threshhold(threshhold)
		{}

		void operator()( VariantIdentifyingData const& snp, Genotypes const& genotypes, Ploidy const& ploidy, genfile::VariantDataReader&, ResultCallback callback ) {
			bool const allDiploid = ( ploidy.array() == 2 ).cast< int >().sum() == ploidy.size() ;
			if( allDiploid ) {
				compute_autosomal_call_mass( snp, genotypes, callback ) ;
			} else {
				return ;
			}
		}
		
		void compute_autosomal_call_mass( VariantIdentifyingData const& snp, Genotypes const& genotypes, ResultCallback callback ) {
			Eigen::VectorXd const masses = genotypes.colwise().sum() ;
			m_hcGenotypes = ( genotypes.array() > m_threshhold ).cast< double >()  * genotypes.array() ;	
			Eigen::VectorXd const hcMasses = m_hcGenotypes.colwise().sum() ; 

			callback( "AA_mass_propn", hcMasses(0)/masses(0) ) ;
			callback( "AB_mass_propn", hcMasses(1)/masses(1) ) ;
			callback( "BB_mass_propn", hcMasses(2)/masses(2) ) ;

			callback( "non-AA-mass_propn", (hcMasses(1)+hcMasses(2))/(masses(1)+masses(2))) ;
			callback( "non-BB-mass_propn", (hcMasses(1)+hcMasses(0))/(masses(1)+masses(0))) ;
		}
		
		std::string get_summary( std::string const& prefix, std::size_t column_width ) const {
			return prefix + "CallMassComputation" ;
		}
	private:
		double const m_threshhold ;
		Genotypes m_hcGenotypes ;
	} ;
	
	struct MissingnessComputation: public SNPSummaryComputation {
		MissingnessComputation( double call_threshhold = 0.9 ): m_call_threshhold( call_threshhold ) {}
		void operator()(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy,
			genfile::VariantDataReader&,
			ResultCallback callback
		) {
			assert( std::size_t( genotypes.rows() ) == ploidy.size() ) ;

			double missingness = double( genotypes.rows() ) - genotypes.array().sum() ;
			callback( "missing_proportion", missingness / double( genotypes.rows() ) ) ;

			bool const allDiploid = ( ploidy.array() == 2 ).cast< int >().sum() == ploidy.size() ;
			if( allDiploid ) {
				callback( "A", 0 ) ;
				callback( "B", 0 ) ;
				if( snp.number_of_alleles() == 2 ) {
					callback( "AA", genotypes.col(0).sum() ) ;
					callback( "AB", genotypes.col(1).sum() ) ;
					callback( "BB", genotypes.col(2).sum() ) ;
				}
				callback( "NULL", genotypes.rows() - genotypes.sum() ) ;
			} else {
				compute_haploid_diploid_counts( snp, genotypes, ploidy, callback ) ;
			}
			callback( "total", genfile::VariantEntry::Integer( genotypes.rows() )) ;
		}
		
		void compute_haploid_diploid_counts(
			VariantIdentifyingData const& snp,
			Genotypes const& genotypes,
			Ploidy const& ploidy,
			ResultCallback callback
		) {
			std::map< int, Eigen::VectorXd > counts ;
			std::map< int, double > null_counts ;
			counts[ -1 ] = Eigen::VectorXd::Zero( 3 ) ;
			counts[ 0 ] = Eigen::VectorXd::Zero( 3 ) ;
			counts[ 1 ] = Eigen::VectorXd::Zero( 3 ) ;
			counts[ 2 ] = Eigen::VectorXd::Zero( 3 ) ;
			std::map< int, std::size_t > sample_counts ;

			for( std::size_t i = 0; i < ploidy.size(); ++i ) {
				counts[ ploidy(i) ] += genotypes.row( i ) ;
				null_counts[ ploidy(i) ] += ( 1 - genotypes.row(i).sum() ) ;
				++sample_counts[ ploidy(i) ] ;
#if DEBUG_SNP_SUMMARY_COMPUTATION
				if( ploidy(i) == 1 && genotypes(i,2) != 0 ) {
					std::cerr << "! ( MissingnessComputation::compute_sex_chromosome_counts() ): individual " << (i+1) << "is male but has genotype " << genotypes.row(i) << "!!\n" ;
				}
#endif
			}
			
			callback( "A", counts[ 1 ]( 0 ) ) ;
			callback( "B", counts[ 1 ]( 1 ) ) ;
			callback( "AA", counts[ 2 ]( 0 ) ) ;
			callback( "AB", counts[ 2 ]( 1 ) ) ;
			callback( "BB", counts[ 2 ]( 2 ) ) ;
			callback( "NULL", null_counts[ 'm' ] + null_counts[ 'f' ] ) ;
			callback( "unknown_ploidy", counts[ -1 ].sum() + null_counts[ -1 ] ) ;
			assert( counts[ 1 ]( 2 ) == 0 ) ;
		}
		
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const { return prefix + "MissingnessComputation" ; }
	private:
		double const m_call_threshhold ;
	} ;

	SNPSummaryComputation::UniquePtr SNPSummaryComputation::create(
		std::string const& name
	) {
		UniquePtr result ;
		if( name == "allele-frequencies" ) { result.reset( new stats::AlleleFrequencyComputation( "everything" )) ; }
		else if( name == "allele-counts" ) { result.reset( new stats::AlleleFrequencyComputation( "counts" )) ; }
		else if( name == "HWE" ) { result.reset( new stats::HWEComputation()) ; }
		else if( name == "missingness" ) { result.reset( new stats::MissingnessComputation()) ; }
		else if( name == "info" ) { result.reset( new stats::InfoComputation()) ; }
		else if( name == "call-mass-proportion" ) { result.reset( new stats::CallMassComputation()) ; }
		else if( name == "intensity-stats" ) { result.reset( new stats::IntensitySummaryComputation() ) ; }
		else if( name == "multi-allele-counts" ) { result.reset( new stats::AlleleCountComputation() ) ; }
		else {
			throw genfile::BadArgumentError( "SNPSummaryComputation::create()", "name=\"" + name + "\"" ) ;
		}
		return result ;
	}
}
