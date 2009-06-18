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

	// Deprecated access functions
	double& proportion_of_AA() { return m_proportion_of_AA ; }
	double& proportion_of_AB() { return m_proportion_of_AB ; }
	double& proportion_of_BB() { return m_proportion_of_BB ; }
	double proportion_of_AA() const { return m_proportion_of_AA ; }
	double proportion_of_AB() const { return m_proportion_of_AB ; }
	double proportion_of_BB() const { return m_proportion_of_BB ; }
	double sum_of_proportions() const { return m_proportion_of_AA + m_proportion_of_BB + m_proportion_of_AB ; }

	void floor() ;
	void round() ;
	
	// New shorter access functions
	double& AA() { return m_proportion_of_AA ; }
	double& AB() { return m_proportion_of_AB ; }
	double& BB() { return m_proportion_of_BB ; }
	double AA() const { return m_proportion_of_AA ; }
	double AB() const { return m_proportion_of_AB ; }
	double BB() const { return m_proportion_of_BB ; }
	double sum() const { return m_proportion_of_AA + m_proportion_of_BB + m_proportion_of_AB ; }

	bool operator==( GenotypeProportions const& other ) const ;

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

