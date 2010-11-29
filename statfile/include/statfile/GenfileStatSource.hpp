#ifndef STATFILE_GENFILESTATSOURCE_HPP
#define STATFILE_GENFILESTATSOURCE_HPP

#include "string_utils.hpp"
#include "genfile/SNPDataSource.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"

namespace statfile {
	// A StatSource which reads from a genfile::SNPDataSource
	class GenfileStatSource: public ColumnNamingStatSource< BuiltInTypeStatSource >, public IstreamAggregator
	{
		typedef ColumnNamingStatSource< BuiltInTypeStatSource > base_t ;
	public:
		GenfileStatSource( std::string const& filename  ):
			m_source( genfile::SNPDataSource::create( genfile::wildcard::find_files_by_chromosome( filename )))
		{
			setup( *m_source )
		}

	public:
		std::size_t number_of_rows() const { return m_source->total_number_of_snps() ; }

		void reset_to_start() {
			m_source->reset_to_start() ;
		}

	protected:

		void read_value( int32_t& value ) {
			if( current_column() == 0 ) {
				get_row() ;
			}
		}

		void read_value( uint32_t& ) {
			if( current_column() == 0 ) {
				get_row() ;
			}
		}
		
		void read_value( std::string& ) {
			if( current_column() == 0 ) {
				get_row() ;
			}
			
		}
		
		void read_value( double& ) {
			if( current_column() == 0 ) {
				get_row() ;
			}
			
		}
		
		void ignore_value() {
			if( current_column() == 0 ) {
			get_row() ;
			}
		}
		
		void ignore_all() {
			if( current_column() == 0 ) {
				get_row() ;
			}
		}

	private:
		genfile::SNPDataSource::UniquePtr m_source ;
		GenRow m_row ;

		void setup( genfile::SNPDataSource const& source ) {
			add_column( "SNPID" ) ;
			add_column( "RSID" ) ;
			add_column( "chromosome" ) ;
			add_column( "position" ) ;
			add_column( "allele1" ) ;
			add_column( "allele2" ) ;
			for( std::size_t i = 0; i < source->number_of_samples(); ++i ) {
				add_column( "AA." + to_string( i ) ) ;
			}
		}
		
		void get_row() {
			m_row.read_from_source( *m_source ) ;
		}
	} ;
	
	
	
}

#endif
