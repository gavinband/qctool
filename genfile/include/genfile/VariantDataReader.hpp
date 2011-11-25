#ifndef GENFILE_VARIANT_DATA_READER_HPP
#define GENFILE_VARIANT_DATA_READER_HPP

#include <memory>
#include <vector>
#include <string>
#include <set>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
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
		typedef genfile::VariantEntry Entry ;
	public:
		struct PerVariantSetter: public vcf::EntriesSetter, public boost::noncopyable {
			typedef std::auto_ptr< PerVariantSetter > UniquePtr ;
			virtual ~PerVariantSetter() throw() {}
			typedef Eigen::MatrixXd Matrix ;
			virtual void operator()( std::auto_ptr< Matrix > value ) ;
		} ;
		struct PerSampleSetter: public vcf::EntriesSetter, public boost::noncopyable {
			typedef std::auto_ptr< PerSampleSetter > UniquePtr ;
			virtual ~PerSampleSetter() throw() {}
			virtual void set_number_of_samples( std::size_t n ) = 0 ;
			virtual void set_sample( std::size_t i ) = 0 ;
		} ;
		typedef boost::function< void ( std::string, std::string ) > SpecSetter ;
	public:
		virtual ~VariantDataReader() {} ;
		virtual VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) = 0 ;
		virtual VariantDataReader& get( std::string const& spec, PerVariantSetter& setter ) ;
		virtual bool supports( std::string const& spec ) const = 0 ;
		virtual void get_supported_specs( SpecSetter ) const = 0 ;
		
		// Convenience method.
		VariantDataReader& get( std::string const& spec, std::vector< std::vector< Entry > >& data ) ;
		// Convenience method setting SingleSNPGenotypeProbabilities.
		VariantDataReader& get( std::string const& spec, SingleSNPGenotypeProbabilities& data ) ;
	} ;
}

#endif
