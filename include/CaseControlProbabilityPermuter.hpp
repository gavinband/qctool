#ifndef CASE_CONTROL_PROBABILITY_PERMUTER_HPP
#define CASE_CONTROL_PROBABILITY_PERMUTER_HPP

#include <vector>
#include "GenotypeProportions.hpp"
#include "CaseControlPermutationsFileReader.hpp"
#include "Timer.hpp"

class CaseControlProbabilityPermuter
{
public:
	CaseControlProbabilityPermuter( std::string const& permutations_filename )
	: 	m_current_permutation(0),
		m_control_genotype_amounts( 0.0, 0.0, 0.0 ),
		m_case_genotype_amounts( 0.0, 0.0, 0.0 ),
		m_genotype_probabilities(0)
	{
		CaseControlPermutationsFileReader file_reader( permutations_filename ) ;
		m_permutations = file_reader.permutations() ;
		m_number_of_controls = file_reader.number_of_zeroes() ;
	}

public:
	std::size_t number_of_permutations() const { return m_permutations.size() ; }
	std::size_t current_permutation() const { return m_current_permutation ; }
	std::size_t size_of_permutations() const { return m_permutations[0].size() ; }
	std::size_t number_of_controls() const { return m_number_of_controls ; }
	std::size_t number_of_cases() const { return size_of_permutations() - m_number_of_controls ; }

	GenotypeAmounts const& get_control_genotype_amounts() const { return m_control_genotype_amounts ; }
	GenotypeAmounts const& get_case_genotype_amounts() const { return m_case_genotype_amounts ; }

	void reset( std::vector< GenotypeProbabilities > const& probabilities ) {
		assert( probabilities.size() == size_of_permutations() ) ;
		m_genotype_probabilities = probabilities ;
		m_current_permutation = 0 ;
		calculate_genotype_amounts_for_permutation( m_current_permutation, &m_control_genotype_amounts, &m_case_genotype_amounts ) ;
	}
	
protected:
	void calculate_genotype_amounts_for_permutation2( std::size_t permutation_i, GenotypeAmounts* control_genotype_amounts, GenotypeAmounts* case_genotype_amounts ) {
		assert( permutation_i < m_permutations.size() ) ;

		GenotypeAmounts control_amounts = GenotypeAmounts( 0.0, 0.0, 0.0 ) ;
		GenotypeAmounts case_amounts = GenotypeAmounts( 0.0, 0.0, 0.0 ) ;

		for( std::size_t i = 0; i < m_permutations[ permutation_i ].size(); ++i ) {
			if( m_permutations[ permutation_i ][i] == 0 ) {
				control_amounts += genotype_probabilities()[i] ;
			}
			else {
				case_amounts += genotype_probabilities()[i] ;
			}
		}
		*control_genotype_amounts = control_amounts ;
		*case_genotype_amounts = case_amounts ;
	}

	void calculate_genotype_amounts_for_permutation( std::size_t permutation_i, GenotypeAmounts* control_genotype_amounts, GenotypeAmounts* case_genotype_amounts ) {
		// Using these stack-allocated GenotypeAmounts is much quicker than using the arguments directly.
		GenotypeAmounts control_amounts = GenotypeAmounts( 0.0, 0.0, 0.0 ) ;
		GenotypeAmounts case_amounts = GenotypeAmounts( 0.0, 0.0, 0.0 ) ;

		// Using these pointers rather than indices seems a lot quicker.
		char const* i = &m_permutations[permutation_i][0] ;
		char const* end_i = &m_permutations[permutation_i][0] + m_permutations[ permutation_i ].size() ;
		GenotypeProbabilities const* genotype_probability_i = &genotype_probabilities()[0] ;
		for( ; i != end_i; ++i, ++genotype_probability_i ) {
			if( *i == 0 ) {
				control_amounts += *genotype_probability_i ;
			}
			else {
				case_amounts += *genotype_probability_i ;
			}
		}
		*control_genotype_amounts = control_amounts ;
		*case_genotype_amounts = case_amounts ;
	}

	std::vector< GenotypeProbabilities > const& genotype_probabilities() const { return m_genotype_probabilities ; }
	std::vector< std::vector< char > > const& permutations() const { return m_permutations ; }

public:
	
	bool move_to_next_permutation() {
		if( (++m_current_permutation) < m_permutations.size() ) {
			calculate_next_genotype_amounts( m_current_permutation, &m_control_genotype_amounts, &m_case_genotype_amounts ) ;
			return true ;
		}
		else {
			m_genotype_probabilities.clear() ;
			return false ;
		}
	}

protected:

	// Given that the passed-in genotype amounts correspond to the amounts for the permutation
	// just before the current one, calculate the amounts for the current permutation.
	virtual void calculate_next_genotype_amounts( std::size_t permutation_i, GenotypeAmounts* control_genotype_amounts, GenotypeAmounts* case_genotype_amounts ) = 0 ;
	
	std::vector< std::vector< char > >& permutations() { return m_permutations ; }

	bool check_first_permutation_has_all_controls_first() const {
		bool found_nonzero = false ;
		for( std::size_t i = 0; i < m_permutations[0].size(); ++i ) {
			found_nonzero = ( m_permutations[ 0 ][ i ] != 0 ) ;
			if( found_nonzero && (i < m_number_of_controls )) {
				return false ;
			}
		}		
		return true ;
	}

private:
	std::vector< std::vector< char > > m_permutations ;
	std::size_t m_size_of_permutations ;
	std::size_t m_current_permutation ;
	std::size_t m_number_of_controls ;
	
	GenotypeAmounts m_control_genotype_amounts, m_case_genotype_amounts ;
	std::vector< GenotypeProbabilities > m_genotype_probabilities ;
} ;


class SimpleCaseControlProbabilityPermuter: public CaseControlProbabilityPermuter
{
public:
	SimpleCaseControlProbabilityPermuter( std::string const& permutations_filename )
		: CaseControlProbabilityPermuter( permutations_filename )
	{}

private:
	void calculate_next_genotype_amounts( std::size_t permutation_i, GenotypeAmounts* control_genotype_amounts, GenotypeAmounts* case_genotype_amounts ) {
		calculate_genotype_amounts_for_permutation( permutation_i, control_genotype_amounts, case_genotype_amounts ) ;
	}
} ;

class DifferentialCaseControlProbabilityPermuter: public CaseControlProbabilityPermuter
{
public:
	DifferentialCaseControlProbabilityPermuter( std::string const& permutations_filename )
		: CaseControlProbabilityPermuter( permutations_filename )
	{
		std::sort( permutations().begin(), permutations().end() ) ;
		if( !check_first_permutation_has_all_controls_first() ) {
			throw PermutationFileFirstPermutationMalformedError() ;
		}
		m_differential_permutations = calculate_differential_permutations() ;
	}

private:

	void calculate_next_genotype_amounts( std::size_t permutation_i, GenotypeAmounts* control_genotype_amounts, GenotypeAmounts* case_genotype_amounts ) {
		GenotypeAmounts case_genotype_amount_adjustment( 0.0, 0.0, 0.0 ) ;

		// Using these pointers rather than indices seems a lot quicker.
		char const* i = &m_differential_permutations[ permutation_i ][0] ;
		char const* end_i = &m_differential_permutations[ permutation_i ][0] + m_differential_permutations[ permutation_i ].size() ;
		GenotypeProbabilities const* genotype_probability_i = &genotype_probabilities()[0] ;
		for( ; i != end_i; ++i, ++genotype_probability_i ) {
			if( *i < 0 ) {
				case_genotype_amount_adjustment -= *genotype_probability_i ;
			}
			else if( *i > 0 ){
				case_genotype_amount_adjustment += *genotype_probability_i ;
			}
		}
		(*control_genotype_amounts) -= case_genotype_amount_adjustment ;
		(*case_genotype_amounts) += case_genotype_amount_adjustment ;
	}

private:

	std::vector< std::vector< char > > calculate_differential_permutations() {
		std::vector< std::vector< char > > differential_permutations( permutations().size() ) ;
		differential_permutations[0] = permutations()[0] ;
		for( std::size_t i = 1; i < number_of_permutations(); ++i ) {
			differential_permutations[i] = calculate_entry_diffences( permutations()[i-1], permutations()[i] ) ;
		}
		return differential_permutations ;
	}

	// Return a vector whose ith entry is vector2[i] - vector1[i].
	// To spell it out, the entry is
	//   0 if the two entries are equal
	//  -1 if vector1[i] is 1 and vector2[i] is 0
	//   1 if vector1[i] is 0 and vector2[i] is 1
	std::vector< char > calculate_entry_diffences( std::vector<char> const& vector1, std::vector<char> const& vector2 ) const {
		assert( vector1.size() == vector2.size() ) ;
		std::vector< char > result( vector1.size() ) ;
		for( std::size_t i = 0; i < vector1.size(); ++i ) {
			result[i] = vector2[i] - vector1[i] ;
		}
		return result ;
	}

private:

	std::vector< std::vector< char > > m_differential_permutations ;
} ;


#endif
