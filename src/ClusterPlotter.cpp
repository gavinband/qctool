
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/string_utils/substitute.hpp"
#include "worker/Worker.hpp"
#include "Eigen/Core"
#include "ClusterPlotter.hpp"
#include "../config.hpp"
#if HAVE_MGL
	#include "mgl/mgl.h"
	#include "mgl/mgl_zb.h"
#endif

void ClusterPlotter::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Cluster plot options" ) ;
#if HAVE_MGL
	options[ "-cluster-plot" ]
		.set_description( "Create a plot of the intensities and genotypes for each SNP. "
		 	"Currently this uses the MathGL library (http://mathgl.sourceforge.net) and so only works "
			"if this is detected during compilation.")
		.set_takes_single_value() ;
	options[ "-cluster-plot-filename" ]
		.set_description( "Template for filename of cluster plots to create.  The following symbols will be replaced with appropriate strings, one per SNP:"
		 	" #analysis (analysis name, set using the -analysis-name option),"
			" #rsid (rsid of SNP),"
			" #position (chromosome and position of SNP),"
			" #intensity (name of intensity field)."
		)
		.set_takes_single_value()
		.set_default_value( "#analysis_#rsid_#position_#intensity.png" ) ;
	options[ "-cluster-plot-intensity-threshhold" ]
		.set_description( "Set the threshhold for intensities in the cluster plots." )
		.set_takes_single_value()
		.set_default_value( 5.0 ) ;
#endif
}

ClusterPlotter::UniquePtr ClusterPlotter::create( appcontext::OptionProcessor const& options, worker::Worker* worker ) {
	UniquePtr result(
		new ClusterPlotter(
			options.get< std::string >( "-analysis-name" ),
			options.get< std::string >( "-cluster-plot-filename" ),
			genfile::string_utils::split_and_strip_discarding_empty_entries( options.get< std::string >( "-cluster-plot" ), ",", " \t" ),
			worker
		)
	) ;
	result->set_intensity_threshhold( options.get< double >( "-cluster-plot-intensity-threshhold" ) ) ;
	return result ;
}

ClusterPlotter::ClusterPlotter(
	std::string const& analysis_name,
	std::string const& filename_template,
	std::vector< std::string > const& call_fields,
	worker::Worker* worker
):
	m_analysis_name( analysis_name ),
	m_filename_template( filename_template ),
	m_call_fields( call_fields ),
	m_intensity_field( "XY" ),
	m_intensity_threshhold( 5.0 ),
	m_worker( worker ),
	m_max_tasks( 10 )
{}

void ClusterPlotter::set_intensity_threshhold( double value ) {
	m_intensity_threshhold = value ;
}

void ClusterPlotter::begin_processing_snps( std::size_t number_of_samples ) {
	m_number_of_samples = number_of_samples ;
}

#if HAVE_MGL

namespace impl {
	struct PlotTask: public worker::Task {
		typedef Eigen::MatrixXd IntensityMatrix ;
		typedef std::vector< int > Genotypes ;
		typedef std::map< std::string, Genotypes > Calls ;
		PlotTask(
			std::string const& analysis_name,
			std::string const& intensity_field,
			std::vector< std::string > const& call_fields,
			genfile::SNPIdentifyingData const& snp,
			genfile::VariantDataReader& data_reader,
			std::string const& filename
		):
			m_intensity_field( intensity_field ),
			m_call_threshhold( 0.9 ),
			m_snp( snp ),
			m_analysis_name( analysis_name ),
			m_filename( filename ),
			m_intensity_threshhold( 5.0 )
		{
			for( std::size_t i = 0; i < call_fields.size(); ++i ) {
				m_calls[ call_fields[i] ] = Genotypes() ;
			}
			genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities ) ;
			data_reader.get( m_intensity_field, intensity_setter ) ;

			for( Calls::iterator i = m_calls.begin(), end_i = m_calls.end(); i != end_i; ++i ) {
				genfile::vcf::ThreshholdingGenotypeSetter< Genotypes > genotype_setter( i->second, m_call_threshhold, 3, 0, 1, 2 ) ;
				data_reader.get( i->first, genotype_setter ) ;
				assert( i->second.size() == std::size_t( m_intensities.rows() ) ) ;
			}
		}

		void set_intensity_threshhold( double const value ) {
			assert( value > 0 ) ;
			m_intensity_threshhold = value ;
		}
		
		void operator()() {
			int M = std::ceil( std::sqrt( m_calls.size() ) ) ;
			int N = std::ceil( m_calls.size() / M ) ;
			mglGraphZB graph( 400 * N, 400 * M ) ;
			graph.Light( true ) ;
			graph.Clf( mglColor( 1.0, 1.0, 1.0 ) ) ;
			mglData x( m_intensities.rows() ), y( m_intensities.rows() ), colour( m_intensities.rows() ) ;
			
			std::size_t count = 0 ;
			for( Calls::iterator i = m_calls.begin(), end_i = m_calls.end(); i != end_i; ++i, ++count ) {
				x.Set( m_intensities.col(0).data(), m_intensities.rows() ) ;
				y.Set( m_intensities.col(1).data(), m_intensities.rows() ) ;
				colour.Set( i->second ) ;
				graph.Title( ( m_analysis_name + ":" + m_snp.get_rsid() ).c_str(), 0, 4 ) ;
				graph.SubPlot( N, M, count ) ;
				graph.SetTickLen( 0.04 ) ;
				double m_max_value = -1 ;
				for( int sample_i = 0; sample_i < m_intensities.rows(); ++sample_i ) {
					for( int col = 0; col < 2; ++col ) {	
						if( m_intensities(sample_i,col) == m_intensities(sample_i,col) ) {
							m_max_value = std::max( m_max_value, m_intensities(sample_i,col) ) ;
						}
					}
				}
				double x_range_max = std::min(
					m_max_value,
					m_intensity_threshhold
				) ;
				double y_range_max = std::min(
					m_max_value,
					m_intensity_threshhold
				) ;
				graph.SetRanges(
					0.0, x_range_max,
					0.0, y_range_max
				) ;
				graph.SetCut( false ) ;
				graph.Axis( "x", true ) ;
				graph.Axis( "y", true ) ;
				graph.CAxis( 0.0, 3.0 ) ;

				graph.Tens( x, y, colour, "RGBk ." ) ;
				graph.Puts( mglPoint( x_range_max / 2, y_range_max * 3 / 5 ), ( i->first + "/" + m_intensity_field ).c_str() ) ;
			}
			graph.WritePNG( m_filename.c_str(), "", false ) ;
		}
	private:
		Calls m_calls ;
		std::string const m_intensity_field ;
		IntensityMatrix m_intensities ;
		double const m_call_threshhold ;
		genfile::SNPIdentifyingData const m_snp ;
		
		std::string const m_analysis_name ;
		std::string const m_call_field ;
		std::string const m_filename ;
		double m_intensity_threshhold ;
	} ;
}

void ClusterPlotter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	using genfile::string_utils::substitute ;
	using genfile::string_utils::to_string ;
	std::string filename = m_filename_template ;
	filename = substitute( filename , "#analysis", m_analysis_name ) ;
	filename = substitute( filename , "#rsid", snp.get_rsid() ) ;
	filename = substitute( filename, "#intensity", m_intensity_field ) ;
	filename = substitute( filename, "#position", to_string( snp.get_position() ) ) ;

	std::auto_ptr< impl::PlotTask > plot_task(
		new impl::PlotTask(
			m_analysis_name,
			m_intensity_field,
			m_call_fields,
			snp,
			data_reader,
			filename
		)
	) ;

	plot_task->set_intensity_threshhold( m_intensity_threshhold ) ;
	
	if( m_tasks.size() < m_max_tasks ) {
		m_tasks.push_back( plot_task ) ;
		m_worker->tell_to_perform_task( m_tasks.back() ) ;
	} else {
		// wait until a task is free...
		std::size_t i = m_tasks.size() ;
		while( i == m_tasks.size() ) {
			for( i = 0; i < m_tasks.size(); ++i ) {
				if( m_tasks[i].check_if_complete() ) {
					break ;
				}
				boost::this_thread::sleep( boost::posix_time::milliseconds( 10 ) ) ;
			}
		}
		assert( i != m_tasks.size() ) ;
		m_tasks.replace( i, plot_task ) ;
		m_worker->tell_to_perform_task( m_tasks[i] ) ;
	}
}
#else
void ClusterPlotter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	assert(0) ;
}
#endif

void ClusterPlotter::end_processing_snps() {
	// nowt to do.
}
