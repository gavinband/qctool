
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_GP_THRESHOLDING_GT_SETTER
#define GENFILE_GP_THRESHOLDING_GT_SETTER

#include "genfile/VariantDataReader.hpp"
#include "genfile/Error.hpp"
#include "genfile/types.hpp"

namespace genfile {
	struct GPThresholdingGTSetter: public VariantDataReader::PerSampleSetter {
		GPThresholdingGTSetter( VariantDataReader::PerSampleSetter& target, double const threshhold ) ;
		~GPThresholdingGTSetter() throw() ;
		
		void initialise( std::size_t nSamples, std::size_t nAlleles ) ;
		bool set_sample( std::size_t i ) ;
		void set_number_of_entries( uint32_t ploidy, std::size_t number_of_entries, OrderType const order_type, ValueType const value_type ) ;
		void set_value( std::size_t, MissingValue const value ) ;
		void set_value( std::size_t, std::string& value ) ;
		void set_value( std::size_t, Integer const value ) ;
		void set_value( std::size_t, double const value ) ;
		void finalise() ;
	private:
		VariantDataReader::PerSampleSetter& m_target ;
		double const m_threshhold ;
		std::size_t m_number_of_alleles ;
		OrderType m_order_type ;
		uint32_t m_ploidy ;
		int m_entry_i ;
		Eigen::VectorXd m_calls ;
		bool m_missing ;
	
	private:
		void send_results() ;
	} ;
}

#endif
