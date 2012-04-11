
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef TrivialSNPDataSink_HPP
#define TrivialSNPDataSink_HPP

#include <iostream>
#include <string>
#include "SNPDataSink.hpp"

namespace genfile {
	// class TrivialSNPDataSink represents a SNPDataSink
	// which never fails and always discards its output.
	class TrivialSNPDataSink: public SNPDataSink
	{
	public:
		operator bool() const {
			return true ;
		}
		
		std::string get_spec() const {
			return "TrivialSNPDataSink" ;
		}

	private:
		void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability,
			Info const& info
		) {
			return ;
		}
	} ;
}

#endif