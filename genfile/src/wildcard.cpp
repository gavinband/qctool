#include <iostream>
#include <string>
#include <map>
#include <cassert>
#include "../config.hpp"
#if HAVE_BOOST_FILESYSTEM
	#include <boost/filesystem.hpp>
	namespace BFS = boost::filesystem ;
#endif
#include "genfile/wildcard.hpp"
#include "genfile/Error.hpp"
#include "genfile/Chromosome.hpp"

namespace genfile {
	namespace wildcard {
		namespace impl {
		
			bool string_has_prefix_and_suffix( std::string const& string_to_check, std::string const& prefix, std::string const& suffix ) {
				if( string_to_check.size() < ( prefix.size() + suffix.size() )) {
					return false ;
				}
				if( prefix != string_to_check.substr( 0, prefix.size()) ) {
					return false ;
				}
				if( suffix != string_to_check.substr( string_to_check.size() - suffix.size(), suffix.size() )) {
					return false ;
				}
				return true ;
			}

#if HAVE_BOOST_FILESYSTEM
			std::vector< FilenameMatch > find_files_with_prefix_and_suffix( BFS::path const& directory, std::string const& prefix, std::string const& suffix ) {
				std::vector< FilenameMatch > result ;
				BFS::directory_iterator dir_i( directory ), end_dir_i ;
				for( ; dir_i != end_dir_i; ++dir_i ) {
#if BOOST_FILESYSTEM_VERSION < 3
					std::string filename = dir_i->filename() ;
#else
					std::string filename = dir_i->path().filename() ;
#endif
					if( string_has_prefix_and_suffix( filename, prefix, suffix )) {
						std::string matching_part = filename.substr( prefix.size(), filename.size() - suffix.size() - prefix.size() ) ;
						result.push_back( FilenameMatch( (directory / filename).string(), matching_part )) ;
					}
				}
				return result ;
			}

			BFS::path get_dir_part( std::string const& path ) {
				BFS::path dir = BFS::path( path ).parent_path() ;
				if( dir.empty() ) {
					dir = "." ;
				}
				return dir ;
			}
#endif
			bool has_wildcard( std::string const& a_string, char const wildcard ) {
				return a_string.find( wildcard ) != std::string::npos ;
			}
		
#if HAVE_BOOST_FILESYSTEM
			std::vector< FilenameMatch > find_files_matching_path_with_chromosomal_wildcard( std::string path, char const wildcard_char ) {
				std::vector< FilenameMatch > result ;
				BFS::path dir = impl::get_dir_part( path ) ;
			  	std::string filename_template = BFS::path( path ).filename() ;
				std::size_t wildcard_pos = filename_template.find( wildcard_char ) ;
				std::string prefix = filename_template.substr( 0, wildcard_pos ) ;
				std::string suffix = filename_template.substr( wildcard_pos + 1, filename_template.size() ) ;

				std::vector< FilenameMatch > candidates = impl::find_files_with_prefix_and_suffix( dir, prefix, suffix ) ;

				// Make a list of the candidates for which the matching part represents a chromosome.
				// Order by the chromosome itself.  Note that the actual matching string may be prefixed with
				// zeroes, so ordering by the matching string is the wrong thing to do.
				std::map< Chromosome, FilenameMatch > matching_files ;
				for( std::size_t i = 0; i < candidates.size(); ++i ) {
					try {
						Chromosome chromosome( candidates[i].match() ) ;
						matching_files[ chromosome ] = candidates[i] ;
					}
					catch( BadArgumentError const& e ) {
						// Filename doesn't indicate a chromosome.
					}
				}
	
				// Copy into the result.
				for( std::map< Chromosome, FilenameMatch >::const_iterator i = matching_files.begin() ;
					i != matching_files.end() ;
					++i
				) {
					result.push_back( i->second ) ;
				}
				return result ;
			}
#endif

			// Look through the filename for any part that appears to specify the chromosome.
			// If we find one, return that chromosome.  Otherwise, return the "Unknown" chromosome
			std::string determine_chromosome_match( std::string const& filename ) {
				for(
					std::size_t pos = filename.find_first_of( "0123456789XYM" ) ;
					pos != std::string::npos && pos != filename.size() ;
					pos = filename.find_first_of( "0123456789XYM", pos )
				) {
					std::size_t pos2 = filename.find_first_not_of( "0123456789XYM", pos ) ;
					if( pos2 == std::string::npos ) {
						pos2 = filename.size() ;
					}
					std::string const possible_match = filename.substr(pos, pos2 - pos) ;
					try {
						Chromosome chromosome( possible_match ) ;
						return chromosome ;
					}
					catch( BadArgumentError const& ) {
						pos = pos2 ;
					}
				}
				// Can't determine chromosome.
				return std::string( Chromosome() ) ;
			}

		}

		bool operator<( FilenameMatch const& left, FilenameMatch const& right ) {
			return ( left.match() < right.match() ) || ( left.filename() < right.filename() ) ;
		}

		bool operator==( FilenameMatch const& left, FilenameMatch const& right ) {
			return ( left.match() == right.match() ) && ( left.filename() == right.filename() ) ;
		}

		std::ostream& operator<<( std::ostream& oStream, FilenameMatch const& match ) {
			return oStream << match.filename() << "(" << match.match() << ")" ;
		}

		namespace impl {
			std::vector< FilenameMatch > choose_chromosomes(
				std::vector< FilenameMatch > matches,
				WildcardMatchChoice choice
			) {
				if( choice != eALL_CHROMOSOMES ) {
					for( std::size_t i = 0; i < matches.size(); ) {
						bool keep_match = false ;
						switch( choice ) {
							case eAUTOSOMAL_CHROMOSOMES:
								keep_match = Chromosome( matches[i].match() ).is_autosome() ;
								break ;
							case eSEX_CHROMOSOMES:
								keep_match = Chromosome( matches[i].match() ).is_sex_determining() ;
								break ;
							case eNON_SEX_CHROMOSOMES:
								keep_match = !Chromosome( matches[i].match() ).is_sex_determining() ;
								break ;
							default:
								assert(0) ;
								break ;
						}
						if( !keep_match ) {
							matches.erase( matches.begin() + i ) ;
						}
						else {
							++i ;
						}
					}
				}
				return matches ;
			}
		}

		std::vector< FilenameMatch > find_files_by_chromosome(
			std::string path,
			WildcardMatchChoice choice,
			char wildcard_char
		) {
			std::vector< FilenameMatch > result ;
#if HAVE_BOOST_FILESYSTEM
			if( BFS::exists( path )) {
				result.push_back( FilenameMatch( path, "" )) ;
				return result ;
			}
			else if( impl::has_wildcard( path, wildcard_char )) {
				result = impl::find_files_matching_path_with_chromosomal_wildcard( path, wildcard_char ) ;
			}
			if( result.empty() ) {
				throw FileNotFoundError( path ) ;
			}

			return impl::choose_chromosomes( result, choice ) ;
#else
			result.push_back( FilenameMatch( path, "" )) ;
			return result ;
#endif
		}
	
		std::vector< FilenameMatch >
		construct_corresponding_filenames(
			std::vector< FilenameMatch > const& input_filename_matches,
			std::string const& output_filename,
			char wildcard_char
		) {
			std::vector< FilenameMatch > result ;

			if( impl::has_wildcard( output_filename, wildcard_char )) {
				std::size_t wildcard_pos = output_filename.find( wildcard_char ) ;
				for( std::size_t i = 0; i < input_filename_matches.size(); ++i ) {
					std::string this_output_filename = output_filename ;
					this_output_filename.replace( wildcard_pos, 1, input_filename_matches[i].match()) ;
					result.push_back( FilenameMatch( this_output_filename, input_filename_matches[i].match())) ;
				}
			}
			else {
				for( std::size_t i = 0; i < input_filename_matches.size(); ++i ) {
					result.push_back( FilenameMatch( output_filename ) ) ;
				}
			}
			return result ;
		}
	}
}
