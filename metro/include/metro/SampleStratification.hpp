
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_SAMPLESTRATIFICATION_HPP
#define CORE_SAMPLESTRATIFICATION_HPP

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "metro/SampleRange.hpp"

namespace metro {
	struct SampleStratification {
		SampleStratification() ;
		SampleStratification( SampleStratification const& other ) ;
		SampleStratification& operator=( SampleStratification const& other ) ;

		void add_stratum( std::string const& name ) ;
		void add_sample( std::string const& strata, int sample_index ) ;
		void add_sample_range( std::string const& strata, metro::SampleRange const& range ) ;

		std::size_t number_of_strata() const { return m_strata_names.size() ; }
		std::size_t size() const { return m_strata_names.size() ; }
		std::string const& stratum_name( std::size_t i ) const ;
		std::vector< metro::SampleRange > const& stratum( std::size_t i ) const ;
		std::vector< metro::SampleRange > const& stratum( std::string const& strata_name ) const ;

	private:
		std::vector< std::string > m_strata_names ;
		std::map< std::string, std::size_t > m_index_by_strata ;
		std::vector< std::vector< metro::SampleRange > > m_sample_ranges ;
	} ;
}

#endif
