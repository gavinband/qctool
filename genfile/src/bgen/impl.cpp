
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#include <stdint.h>
#include <map>
#include "genfile/bgen/impl.hpp"

namespace genfile {
	namespace bgen {
		namespace impl {
			uint32_t get_flags( std::string const& version ) {
				if( version == "v11" ) {
					return e_CompressedSNPBlocks | e_v11Layout ;
				} else if( version == "v12" ) {
					return e_CompressedSNPBlocks | e_v12Layout ;
				} else {
					assert(0) ;
				}
			}
			
			namespace v12 {
				namespace {
					struct CompareFractionalPart{
						CompareFractionalPart( double* v, std::size_t n ):
							m_v( v ), m_n( n )
						{}
						CompareFractionalPart( CompareFractionalPart const& other ):
							m_v( other.m_v ), m_n( other.m_n )
						{}
						CompareFractionalPart& operator=( CompareFractionalPart const& other ) {
							m_v = other.m_v ;
							m_n =  other.m_n ;
							return *this ;
						}
						
						bool operator()( std::size_t a, std::size_t b ) const {
							return ( m_v[a] - std::floor( m_v[a] ) ) > ( m_v[b] - std::floor( m_v[b] )) ;
						}
					private:
						double* m_v ;
						std::size_t m_n ;
					} ;
				}
				void round_probs_to_scaled_simplex( double* p, std::size_t* index, std::size_t const n, int const number_of_bits ) {
					double const scale = ( 0xFFFFFFFFFFFFFFFF >> ( 64 - number_of_bits ) ) ;
					std::multimap< double, double*, std::greater< double > > fractional_parts ;
					double total_fractional_part = 0.0 ;
					for( std::size_t i = 0; i < n; ++i ) {
						p[i] *= scale ;
						index[i] = i ;
						total_fractional_part += p[i] - std::floor( p[i] ) ;
					}
					std::size_t const upper = std::floor( total_fractional_part + 0.5 ) ;
					std::sort( index, index + n, CompareFractionalPart( p, n ) ) ;

	#if DEBUG_BGEN_FORMAT > 2
		std::cerr << "round_probs_to_scaled_simplex(): number_of_bits = " << number_of_bits
				<< ", scale = " << scale
				<< ", total_fractional_part = " << total_fractional_part
				<< ", upper = " << upper << ".\n" ;
		std::cerr << "round_probs_to_scaled_simplex(): p1 = " << *p << ".\n" ;
	#endif

					for( std::size_t i = 0; i < n; ++i ) {
						if( i < upper ) {
							p[ index[i] ] = std::ceil( p[ index[i] ] ) ;
						} else {
							p[ index[i] ] = std::floor( p[ index[i] ] ) ;
						}
					}
				}

				char* write_scaled_probs(
					uint64_t* data,
					std::size_t* offset,
					double const* probs,
					std::size_t const n,
					int const number_of_bits,
					char* buffer,
					char* const end
				) {
					for( std::size_t i = 0; i < (n-1); ++i ) {
						uint64_t const storedValue = uint64_t( probs[i] ) ;
						*data |= storedValue << (*offset) ;
						(*offset) += number_of_bits ;
						if( (*offset) >= 32 ) {
#if DEBUG_BGEN_FORMAT > 1
							std::cerr << "genfile::bgen:impl::v12::write_scaled_probs(): data = " << std::hex << (*data)
								<< ", buffer = " << reinterpret_cast< void* >( buffer )
								<< std::dec << ", offset = " << (*offset)
								<< ", number_of_bits = " << number_of_bits
								<< ".\n";
							std::cerr << "genfile::bgen:impl::v12::write_scaled_probs(): writing data.\n" ;
#endif
							assert( (buffer+4) <= end ) ;
							buffer = std::copy(
								reinterpret_cast< char const* >( data ),
								reinterpret_cast< char const* >( data ) + 4,
								buffer
							) ;
							(*offset) -= 32 ;
							(*data) >>= 32 ;
#if DEBUG_BGEN_FORMAT > 1
							std::cerr << "genfile::bgen:impl::v12::write_scaled_probs(): after write, data = " << std::hex << (*data)
								<< ", buffer = " << std::hex << reinterpret_cast< void* >( buffer )
								<< ", last four bytes written: " << genfile::string_utils::to_hex( buffer-4, buffer ) << ".\n" ;
#endif
							
						}
					}
					return buffer ;
				}
			}
		}
	}
}
