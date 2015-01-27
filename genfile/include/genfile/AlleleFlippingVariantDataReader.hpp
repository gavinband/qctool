
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_ALLELEFLIPPINGVARIANTDATAREADER_HPP
#define GENFILE_ALLELEFLIPPINGVARIANTDATAREADER_HPP


#include <string>
#include <vector>
#include <algorithm>
#include <boost/format.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	class AlleleFlippingVariantDataReader: public VariantDataReader {
	public:
		AlleleFlippingVariantDataReader(
			std::size_t number_of_samples,
			VariantDataReader::UniquePtr base_reader,
			char flip
		):
			m_number_of_samples( number_of_samples ),
			m_base_reader( base_reader ),
			m_flip( flip )
		{}
	
		std::size_t get_number_of_samples() const { return m_number_of_samples ; }
	
		AlleleFlippingVariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
			if(
				m_flip == StrandAligningSNPDataSource::eNoFlip
			) {
				// fall through to base reader.
				m_base_reader->get( spec, setter ) ;
			}
			else {
				OffsetFlippedAlleleSetter flipped_setter( setter, m_number_of_samples, m_flip, 0 ) ;
				m_base_reader->get( spec, flipped_setter ) ;
			}
			return *this ;
		}
	
		bool supports( std::string const& spec ) const {
			return m_base_reader->supports( spec ) ;
		}
	
		void get_supported_specs( SpecSetter setter ) const {
			return m_base_reader->get_supported_specs( setter ) ;
		}

	private:
		std::size_t const m_number_of_samples ;
		VariantDataReader::UniquePtr m_base_reader ;
		char const m_flip ;
	} ;
}

#endif
	