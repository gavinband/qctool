#ifndef __GTOOL_ALLELEPROPORTIONS__
#define __GTOOL_ALLELEPROPORTIONS__

#include <vector>

/* Example usage:
GenotypeProportions acc( 0.0, 0.0, 0.0 ) ;

acc = std::accumulate( row.begin_allele_proportions(), row.end_allele_proportions(), acc ) ;
*/

// Class 
struct GenotypeProportions
{
    GenotypeProportions() ;
	GenotypeProportions( double aa, double ab, double bb );

	// Modify the contained proportions
	void floor() ;
	void round() ;
	void zero() ;
	
	// New shorter access functions
	double& AA() { return m_proportion_of_AA ; }
	double& AB() { return m_proportion_of_AB ; }
	double& BB() { return m_proportion_of_BB ; }
	double AA() const { return m_proportion_of_AA ; }
	double AB() const { return m_proportion_of_AB ; }
	double BB() const { return m_proportion_of_BB ; }
	double sum() const { return m_proportion_of_AA + m_proportion_of_BB + m_proportion_of_AB ; }

	bool operator==( GenotypeProportions const& other ) const ;
	bool operator!=( GenotypeProportions const& other ) const ;

	// Convenient operators
	GenotypeProportions& operator+=( GenotypeProportions const& right ) {
		m_proportion_of_AA += right.m_proportion_of_AA ;
		m_proportion_of_BB += right.m_proportion_of_BB ;
		m_proportion_of_AB += right.m_proportion_of_AB ;
        return *this ;
	}

	GenotypeProportions& operator/=( double scalar ) {
		m_proportion_of_AA /= scalar ;
		m_proportion_of_BB /= scalar ;
		m_proportion_of_AB /= scalar ;
        return *this ;
	}
	
	private:

		double m_proportion_of_AA ;
		double m_proportion_of_AB ;
		double m_proportion_of_BB ;
} ;


// non-member operators
GenotypeProportions operator+( GenotypeProportions const& left, GenotypeProportions const& right ) ;
GenotypeProportions operator/( GenotypeProportions const& left, double right ) ;
std::ostream& operator<<( std::ostream&, GenotypeProportions const& ) ; 

// Operators on vectors of GenotypeProportions
std::vector< GenotypeProportions > operator+( std::vector<GenotypeProportions> const&, std::vector<GenotypeProportions> const& ) ;
// Divide a vector of GenotypeProportions by a scalar
std::vector< GenotypeProportions > operator/( std::vector<GenotypeProportions> const&, double ) ;

#endif

