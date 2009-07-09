#ifndef __GTOOL_GENROWFILESINK_HPP
#define __GTOOL_GENROWFILESINK_HPP

#include "GenRow.hpp"
#include "SimpleFileObjectSink.hpp"
#include "FileUtil.hpp"
#include "SNPDataSink.hpp"

struct SNPDataSinkGenRowSink: public ObjectSink< GenRow >
{
public:
	SNPDataSinkGenRowSink( std::auto_ptr< genfile::SNPDataSink > snp_data_sink ) ;
	operator bool() const { return *m_snp_data_sink ; }
	SNPDataSinkGenRowSink& write( GenRow const& row ) ;
	
protected:

		std::auto_ptr< genfile::SNPDataSink > m_snp_data_sink ;
	
} ;

#endif
