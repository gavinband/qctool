
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_FILE_OUTPUTTER_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_FILE_OUTPUTTER_HPP

#include <string>
#include <memory>
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"

struct FileOutputter: public boost::noncopyable {
	typedef std::auto_ptr< FileOutputter > UniquePtr ;
	typedef boost::shared_ptr< FileOutputter > SharedPtr ;
	
	static UniquePtr create( std::string const& filename ) { return UniquePtr( new FileOutputter( filename ) ) ; }
	static SharedPtr create_shared( std::string const& filename ) { return SharedPtr( new FileOutputter( filename ) ) ; }

	FileOutputter( std::string const& filename ):
		m_filename( filename ),
		m_sink( statfile::BuiltInTypeStatSink::open( filename ))
	{
		(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "computation_name" | "variable" | "value" ;
	}

	void operator()(
		std::size_t index,
		genfile::SNPIdentifyingData const& snp,
		std::string const& computation_name,
		std::string const& value_name,
		genfile::VariantEntry const& value
	) {
		(*m_sink)
			<< snp.get_SNPID()
			<< snp.get_rsid()
			<< std::string( snp.get_position().chromosome() )
			<< snp.get_position().position()
			<< snp.get_first_allele()
			<< snp.get_second_allele()
			<< computation_name
			<< value_name ;
		if( value == value ) {
			(*m_sink) << value ;
		}
		else {
			(*m_sink) << "NA" ;
		}
		(*m_sink) << statfile::end_row() ;
		;
	}

private:
	std::string const m_filename ;
	statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
} ;


#endif
