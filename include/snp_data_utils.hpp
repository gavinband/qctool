#ifndef SNP_DATA_UTILS_HPP
#define SNP_DATA_UTILS_HPP

namespace gen {
	struct Ignorer
	{
		template< typename T > void operator()( T const& ) {} ;
		template< typename T > void operator()( T const&, T const& ) {} ;
		template< typename T > void operator()( T const&, T const&, T const& ) {} ;
		template< typename T1, typename T2 > void operator()( T1 const&, T2 const&, T2 const&, T2 const& ) {} ;
	} ;

	Ignorer ignore() ;

	template< typename T >
	struct ValueSetter
	{
		ValueSetter( T& t ): m_t(t) {}
		void operator()( T const& t ) { m_t = t ; }
	private:
		T& m_t ;
	} ;

	template< typename T > ValueSetter< T > set_value( T& t ) {
		return ValueSetter< T >( t ) ;
	}
	
}

#endif