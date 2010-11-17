#ifndef GENFILE_CHROMOSOME_HPP
#define GENFILE_CHROMOSOME_HPP

#include <string>
#include <iostream>

namespace genfile {
	
	enum ChromosomeEnum {
		Chromosome1 = 1,
		Chromosome2,
		Chromosome3,
		Chromosome4,
		Chromosome5,
		Chromosome6,
		Chromosome7,
		Chromosome8,
		Chromosome9,
		Chromosome10,
		Chromosome11,
		Chromosome12,
		Chromosome13,
		Chromosome14,
		Chromosome15,
		Chromosome16,
		Chromosome17,
		Chromosome18,
		Chromosome19,
		Chromosome20,
		Chromosome21,
		Chromosome22,
		XChromosome = 23,
		YChromosome = 24,
		XYPseudoAutosomalDNA = 253,
		MitochondrialDNA = 254,
		UnidentifiedChromosome = 255
	} ;
	
	struct Chromosome
	{
		Chromosome()
			: m_chromosome_e( UnidentifiedChromosome )
		{}

		Chromosome( ChromosomeEnum chromosome_e )
			: m_chromosome_e( chromosome_e )
		{}

		Chromosome( unsigned char c )
			: m_chromosome_e( ChromosomeEnum( c ) )
		{}

		Chromosome( std::string const& chromosome_str ) ;
		
		Chromosome& operator=( ChromosomeEnum chromosome_e ) {
			m_chromosome_e = chromosome_e ;
			return *this ;
		}
		
		bool operator==( Chromosome other ) {
			return m_chromosome_e == other.m_chromosome_e ;
		}
		
		bool operator<=( Chromosome other ) {
			return m_chromosome_e <= other.m_chromosome_e ;
		}

		bool operator>=( Chromosome other ) {
			return m_chromosome_e >= other.m_chromosome_e ;
		}

		bool operator<( Chromosome other ) {
			return m_chromosome_e < other.m_chromosome_e ;
		}

		bool operator>( Chromosome other ) {
			return m_chromosome_e > other.m_chromosome_e ;
		}
		
		Chromosome& operator++() ;
		
		bool is_sex_determining() const ;
        
		operator ChromosomeEnum() const { return m_chromosome_e ; }
		
		operator std::string () const ;

	private:
		ChromosomeEnum m_chromosome_e ;
	} ;
	
	std::ostream& operator<<( std::ostream&, Chromosome const& ) ;
	std::istream& operator>>( std::istream&, Chromosome& ) ;


	struct ChromosomeMismatchError: public std::exception
	{
		ChromosomeMismatchError( Chromosome expected, Chromosome got )
			: m_expected( expected ),
			m_got( got )
		{}
		
		virtual ~ChromosomeMismatchError() throw() {}
		char const* what() const throw() { return "ChromosomeMismatchError" ; }

		Chromosome const& expected() const { return m_expected ; }
		Chromosome const& got() const { return m_got ; }
	private:
		Chromosome m_expected, m_got ;
	} ;
}

#endif
