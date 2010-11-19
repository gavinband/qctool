#ifndef APPCONTEXT_OSTREAM_TEE_HPP
#define APPCONTEXT_OSTREAM_TEE_HPP

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <cassert>

namespace appcontext {
	class OstreamTee: public std::ostream  {
	public:
		~OstreamTee() ;

		void add_stream( std::string const& name, std::ostream& stream ) ; 
		void add_stream( std::string const& name, std::auto_ptr< std::ostream > stream ) ; 
		std::ostream& operator[]( std::string const& name ) ;

		typedef std::ostream& (*Manipulator)( std::ostream& ) ;

		template< typename T >
		friend OstreamTee& operator<<( OstreamTee & ostream_tee, T const& t ) ;

		friend OstreamTee& operator<<( OstreamTee & ostream_tee, Manipulator t ) ;

	private:
		std::map< std::string, std::ostream* > m_streams ;
		std::vector< std::ostream* > m_managed_streams ;
	} ;

	template< typename T >
	OstreamTee& operator<<( OstreamTee& ostream_tee, T const& t ) {
		for(
			std::map< std::string, std::ostream* >::const_iterator i = ostream_tee.m_streams.begin() ;
			i != ostream_tee.m_streams.end() ;
			++i
		) {
			(*(i->second)) << t ;
		}

		return ostream_tee ;
	}

}
#endif
