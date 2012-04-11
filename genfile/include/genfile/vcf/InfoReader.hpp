
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_INFO_READER_HPP
#define GENFILE_VCF_INFO_READER_HPP

#include <string>
#include <map>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/function.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/vcf/Types.hpp"

namespace genfile {
	namespace vcf {
		class InfoReader
		//
		// This class reads the INFO fields from a VCF file.
		// Usage:
		//
		// InfoReader( "NS=3;DP=14;AF=0.5;DB;H2", entry_types )
		//		( "NS", set_number_of_samples )
		//		( "DB", set_dbSNP_membership ) ;
		{
		public:
			InfoReader(
				std::size_t number_of_alleles,
				std::string const& data,
				boost::ptr_map< std::string, VCFEntryType > const& entry_types
			) ;
			
			typedef EntriesSetter Setter ;

			InfoReader& operator()( std::string const& spec, Setter& setter ) ;

			// Deconstruct.  This performs the actual work of setting values.
			~InfoReader() ;

		private:
			std::size_t const m_number_of_alleles ;
			std::map< std::string, std::string > const m_data ;
			boost::ptr_map< std::string, VCFEntryType > const& m_entry_types ;
			typedef std::map< std::string, Setter* > Setters ;
			Setters m_setters ;
		private:
			void set_values() const ;
			InfoReader( InfoReader const& other ) ;
			InfoReader& operator=( InfoReader const& other ) ;
		} ;
	}
}

#endif
