#ifndef GENFILE_VCF_GET_SET_EIGEN_HPP
#define GENFILE_VCF_GET_SET_EIGEN_HPP

#include <Eigen/Core>
#include "genfile/vcf/get_set.hpp"

namespace genfile {
	namespace vcf {
		// Genotype setter which stores its genotypes in an Eigen matrix with
		// 3 rows and n columns (where n is the number of samples.)
		template< typename Derived >
		struct GenotypeSetter< Eigen::MatrixBase< Derived > >: public GenotypeSetterBase
		{
			GenotypeSetter( Eigen::MatrixBase< Derived >& result ): m_result( result ) {} ;
			void set_number_of_samples( std::size_t n ) {
				m_result.derived().resize( n, 3 ) ;
				GenotypeSetterBase::set_number_of_samples( n ) ;
			}
			void set( std::size_t sample_i, double AA, double AB, double BB ) {
				m_result( sample_i, 0 ) = AA ;
				m_result( sample_i, 1 ) = AB ;
				m_result( sample_i, 2 ) = BB ;
			}
		private:
			Eigen::MatrixBase< Derived >& m_result ;
		} ;
	}
}

#endif
