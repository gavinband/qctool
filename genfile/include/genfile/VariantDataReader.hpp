
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VARIANT_DATA_READER_HPP
#define GENFILE_VARIANT_DATA_READER_HPP

#include <memory>
#include <vector>
#include <string>
#include <set>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "genfile/VariantEntry.hpp"
#include "genfile/get_set.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/BasicTypes.hpp"
#include "genfile/vcf/Types.hpp"

namespace genfile {
	class VariantDataReader: public boost::noncopyable
	{
	public:
		typedef std::auto_ptr< VariantDataReader > UniquePtr ;
		typedef boost::shared_ptr< VariantDataReader > SharedPtr ;
		typedef genfile::VariantEntry Entry ;
	public:
		struct PerVariantSetter: public vcf::EntriesSetter, public boost::noncopyable {
			typedef std::auto_ptr< PerVariantSetter > UniquePtr ;
			virtual ~PerVariantSetter() throw() {}
			typedef Eigen::MatrixXd Matrix ;
			virtual void operator()( std::auto_ptr< Matrix > value ) ;
		} ;

		typedef vcf::PerSampleEntriesSetter PerSampleSetter ;
		typedef boost::function< void ( std::string, std::string ) > SpecSetter ;
	public:
		virtual ~VariantDataReader() {} ;
		virtual VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) = 0 ;
		// The sole purpose of the next method is to support temporary setters.
		// The method casts away const and forwards to the non-const version.
		// Do NOT use this with truly const setter objects!
		VariantDataReader& get( std::string const& spec, PerSampleSetter const& setter ) ;
		virtual VariantDataReader& get( std::string const& spec, PerVariantSetter& setter ) ;
		virtual bool supports( std::string const& spec ) const = 0 ;
		virtual void get_supported_specs( SpecSetter ) const = 0 ;
		virtual std::size_t get_number_of_samples() const = 0 ;

		// Convenience method.
		VariantDataReader& get( std::string const& spec, std::vector< std::vector< Entry > >& data ) ;
		// Convenience method setting SingleSNPGenotypeProbabilities.
		//VariantDataReader& get( std::string const& spec, SingleSNPGenotypeProbabilities& data ) ;
	} ;
}

#endif
