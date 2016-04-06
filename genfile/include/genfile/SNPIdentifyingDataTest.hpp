
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
	class VariantIdentifyingDataTest
	{
	public:
		typedef std::auto_ptr< VariantIdentifyingDataTest > UniquePtr ;
		typedef boost::shared_ptr< VariantIdentifyingDataTest > SharedPtr ;
	public:
		virtual ~VariantIdentifyingDataTest() {} ;
		
		virtual std::string display() const = 0 ;
	
		// Return true if the variant passes the test.
		virtual bool operator()( VariantIdentifyingData const& data ) const = 0 ;
		
		// Return a vector of indices of SNPs which pass the test.
		std::vector< std::size_t > get_indices_of_filtered_in_snps( std::vector< VariantIdentifyingData> const& snps ) const ;
		typedef boost::function< VariantIdentifyingData const& ( std::size_t ) > SNPGetter ;
		std::vector< std::size_t > get_indices_of_filtered_in_snps( std::size_t number_of_snps, SNPGetter ) const ;
		
	protected:
		VariantIdentifyingDataTest() {} ;
	private:
		// forbid copying, assignment.
		VariantIdentifyingDataTest( VariantIdentifyingDataTest const& ) ;
		VariantIdentifyingDataTest& operator=( VariantIdentifyingDataTest const& ) ;
	} ;
	
	std::ostream& operator<<( std::ostream&, VariantIdentifyingDataTest const& ) ;
	
	struct VariantIdentifyingDataTestNegation: public VariantIdentifyingDataTest
	{
		VariantIdentifyingDataTestNegation( VariantIdentifyingDataTest::UniquePtr ) ;
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	private:
		VariantIdentifyingDataTest::UniquePtr m_subtest ;
	} ;
	
	struct CompoundVariantIdentifyingDataTest: public VariantIdentifyingDataTest
	{
		virtual ~CompoundVariantIdentifyingDataTest() ;
		void add_subtest( VariantIdentifyingDataTest::UniquePtr subtest ) ;
		std::size_t get_number_of_subtests() const ;
		VariantIdentifyingDataTest const& get_subtest( std::size_t index ) const ;
	private:
		std::vector< VariantIdentifyingDataTest* > m_subtests ;
	} ;

	struct VariantIdentifyingDataTestConjunction: public CompoundVariantIdentifyingDataTest
	{
	public:
		typedef std::auto_ptr< VariantIdentifyingDataTestConjunction > UniquePtr ;
	public:
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	} ;

	struct VariantIdentifyingDataTestDisjunction: public CompoundVariantIdentifyingDataTest
	{
	public:
		typedef std::auto_ptr< VariantIdentifyingDataTestDisjunction > UniquePtr ;
	public:
		bool operator()( VariantIdentifyingData const& data ) const ;
		std::string display() const ;
	} ;
}

#endif
