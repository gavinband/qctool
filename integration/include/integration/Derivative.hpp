#ifndef INTEGRATION_DERIVATIVE_HPP
#define INTEGRATION_DERIVATIVE_HPP

namespace integration {
	
	template< typename Function >
	struct Derivative
	{
		typedef typename Function::Vector Vector ;
		typedef typename Function::Matrix Matrix ;
		
		Derivative( Function& function ): m_function( function ) {}
		Derivative( Derivative const& other ): m_function( other.m_function ) {}
		
		void evaluate_at( Vector const& parameters ) { m_function.evaluate_at( parameters ) ; }
		Vector get_value_of_function() const { return m_function.get_value_of_first_derivative() ; }
		Matrix get_value_of_first_derivative() const { return m_function.get_value_of_second_derivative() ; }
	private:
		Function& m_function ;
		Derivative& operator=( Derivative const& other ) ;
	} ;
	
	template< typename Function >
	Derivative< Function > derivative( Function& function ) {
		return Derivative< Function >( function ) ;
	}
}

#endif
