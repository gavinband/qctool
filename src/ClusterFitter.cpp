#include "genfile/SNPDataSourceProcessor.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "ClusterFitter.hpp"
#include "NormalClusterFitter.hpp"

void ClusterFitter::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Cluster options" ) ;
	options[ "-fit-clusters" ]
		.set_description( "Using intensity and genotype data, re-fit clusters for each SNP and output them. "
		 	"The argument specifies the key of the genotype and intensity fields to work with. "
			"This should be a comma-separated list of values of the form \"genotypes/intensities\"." )
		.set_takes_single_value()
		.set_default_value( "genotypes/intensities" ) ;
		
	options[ "-fit-cluster-file" ]
		.set_description( "Override the default output file for cluster fits.")
		.set_takes_single_value()
		.set_default_value( "qctool.clusters" ) ;
	options[ "-call-threshhold" ]
		.set_description( "Set the threshhold for making genotype calls from call probabilities, where applicable." )
		.set_takes_single_value()
		.set_default_value( 0.9 ) ;
}

ClusterFitter::UniquePtr ClusterFitter::create( appcontext::OptionProcessor const& options, genfile::SNPDataSourceProcessor::Callback* output ) {
	assert( !output ) ;
	return UniquePtr( new NormalClusterFitter( options ) ) ;
}

