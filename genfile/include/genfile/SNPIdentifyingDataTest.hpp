
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNPIDENTIFYINGDATATEST_HPP
#define GENFILE_SNPIDENTIFYINGDATATEST_HPP

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "unistd.h"
#include "genfile/GenomePosition.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantIdentifyingData.hpp"

namespace genfile {
	
	// Base class for objects which represent a boolean test of a SNP's "identifying data"
	// (i.e. position, ids, and alleles)
	//
	class SNPIdentifyingDataTest
	{
	public:
		typedef std::auto_ptr< SNPIdentifyingDataTest > UniquePtr ;
		typedef boost::shared_ptr< SNPIdentifyingDataTest > SharedPtr ;
	public:
		virtual ~SNPIdentifyingDataTest() {} ;
		
		virtual std::string display() const = 0 ;
	
		// Return true if the variant passes the test.
		virtual bool operator()( VariantIdentifyingData const& data ) const = 0 ;
		
		// Return a vector of indices of SNPs which pass the test.
		std::vector< std::size_t > get_indices_of_filtered_in_snps( std::vector< VariantIdentifyingData> const& snps ) const ;
		typedef boost::function< VariantIdentifyingData const& ( std::size_t ) > SNPGetter ;
		std::vector< std::size_t > get_indices_of_filtered_in_snps( std::size_t number_of_snps, SNPGetter ) const ;
		
	protected:
		SNPIdentifyingDataTest() {} ;
	private:
		// forbid copying, assignment.
		SNPIdentifyingDataTest( SNPIdentifyingDataTest const& ) ;
		SNPIdentifyingDataTest& operator=( SNPIdentifyingDataTest const& ) ;
	} ;
	
	std::ostream& operator<<( std::ostream&, SNPIdentifyingDataTest const& ) ;
	
	struct SNPIdentifyingDataTestNegation: public SNPIdentifyingDataTest
	{
		SNPIdentifyingDataTestNegation( SNPIdentifyingDataTest::UniquePtr ) ;
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	private:
		SNPIdentifyingDataTest::UniquePtr m_subtest ;
	} ;
	
	struct CompoundSNPIdentifyingDataTest: public SNPIdentifyingDataTest
	{
		virtual ~CompoundSNPIdentifyingDataTest() ;
		void add_subtest( SNPIdentifyingDataTest::UniquePtr subtest ) ;
		std::size_t get_number_of_subtests() const ;
		SNPIdentifyingDataTest const& get_subtest( std::size_t index ) const ;
	private:
		std::vector< SNPIdentifyingDataTest* > m_subtests ;
	} ;

	struct SNPIdentifyingDataTestConjunction: public CompoundSNPIdentifyingDataTest
	{
	public:
		typedef std::auto_ptr< SNPIdentifyingDataTestConjunction > UniquePtr ;
	public:
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	} ;

	struct SNPIdentifyingDataTestDisjunction: public CompoundSNPIdentifyingDataTest
	{
	public:
		typedef std::auto_ptr< SNPIdentifyingDataTestDisjunction > UniquePtr ;
	public:
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	} ;
}

#endif
