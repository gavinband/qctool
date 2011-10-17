#include <cassert>
#include <vector>
#include <Eigen/Core>
#include "NormalClusterIntensityModel.hpp"
#include "MultiVariateNormalDistribution.hpp"

NormalClusterIntensityModel::NormalClusterIntensityModel() {
}

void NormalClusterIntensityModel::add_cluster( Vector const& mean, Matrix const& variance, double weight ) {
	std::auto_ptr< Model::ComponentDistribution > cluster(
		new MultivariateNormalDistribution< Vector, Matrix >( mean, variance )
	) ;
	m_model.add_component( cluster, weight ) ;
}

std::size_t NormalClusterIntensityModel::get_number_of_genotype_classes() const {
	return m_model.get_number_of_components() ;
}

double NormalClusterIntensityModel::get_log_density_at_intensity( Vector const& v ) const {
	return m_model.get_log_density( v ) ;
}

double NormalClusterIntensityModel::get_log_density_at_intensity_given_genotype( std::size_t genotype, Vector const& v ) const {
	return m_model.get_log_density_given_component( genotype, v ) ;
}

