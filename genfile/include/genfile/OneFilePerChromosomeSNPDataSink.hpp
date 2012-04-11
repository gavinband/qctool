
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_ONEFILEPERCHROMOSOME_SNP_DATA_SINK_HPP
#define GENFILE_ONEFILEPERCHROMOSOME_SNP_DATA_SINK_HPP

#include <memory>
#include <map>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPDataSink.hpp"

namespace genfile {
	class OneFilePerChromosomeSNPDataSink: public SNPDataSink
	{
	public:
		typedef std::auto_ptr< OneFilePerChromosomeSNPDataSink > UniquePtr ;
		static UniquePtr create(
			std::string const& filename,
			Metadata const& metadata = Metadata(),
			std::string const& wildcard = "#"
		) ;
		
		OneFilePerChromosomeSNPDataSink(
			std::string const& filename,
			Metadata const& metadata = Metadata(),
			std::string const& wildcard = "#"
		) ;

	public:
		operator bool() const ;
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
		) ;

		std::string get_spec() const ;

	private:
		boost::tuple< std::string, std::string, std::string > m_filename_template ;
		Metadata const m_metadata ;
		std::map< Chromosome, SNPDataSink* > m_sinks ;
		boost::ptr_vector< SNPDataSink > m_sink_storage ;
		std::map< Chromosome, std::string > m_filenames ;
		
		SNPDataSink& get_sink_for_chromosome( Chromosome const& chromosome ) ;
	} ;
}

#endif
