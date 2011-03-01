#include <cmath>
#include <utility>
#include "integration/Integrator.hpp"
#include "integration/GanderGautschiAdaptiveIntegratorBase.hpp"

namespace integration {

	GanderGautschiAdaptiveIntegratorBase::GanderGautschiAdaptiveIntegratorBase( double desired_error ): Integrator( desired_error )
	{
		// A maximum of 13 function evaluations, divided into symmetric pairs are used by the three quadrature rules.
		// Only the even-numbered of these points are used except for a once-time global integral estimate.
		m_evaluation_points.resize( 7 ) ;
		m_evaluation_points[0] = 1 ;
		m_evaluation_points[1] = .94288241569547971905635175843185720232 ; // x1
		m_evaluation_points[2] = std::sqrt( 2.0 / 3.0 ) ;
		m_evaluation_points[3] = .64185334234578130578123554132903188354 ; // x2
		m_evaluation_points[4] = 1.0 / std::sqrt( 5.0 ) ;
		m_evaluation_points[5] = .23638319966214988028222377349205292599 ; // x3
		m_evaluation_points[6] = 0.0 ;
		
		// We list the quadrature rules as (point, weight).
		m_rules.resize( 3 ) ;
		m_rules[0].resize( 2 ) ;
		m_rules[0][0] = std::make_pair( 0,	1.0 / 6.0 ) ;
		m_rules[0][1] = std::make_pair( 4, 	5.0 / 6.0 ) ;

		m_rules[1].resize( 4 ) ;
		m_rules[1][0] = std::make_pair( 0, 	11.0 / 210.0 ) ;
		m_rules[1][1] = std::make_pair( 2, 	72.0 / 245.0 ) ;
		m_rules[1][2] = std::make_pair( 4, 	125.0 / 294.0 ) ;
		m_rules[1][3] = std::make_pair( 6, 	16.0 / 35.0 ) ;

		m_rules[2].resize( 7 ) ;
		m_rules[2][0] = std::make_pair( 0,	.015827191973480183087169986733305510591 ) ;
		m_rules[2][1] = std::make_pair( 1,	.094273840128850045531282505077108171960 ) ;
		m_rules[2][2] = std::make_pair( 2,	.15507198733658539625363597980210298680 ) ;
		m_rules[2][3] = std::make_pair( 3,	.18882157396018245442000533937297167125 ) ;
		m_rules[2][4] = std::make_pair( 4,	.19977340522685852679206802206648840246 ) ;
		m_rules[2][5] = std::make_pair( 5,	.22492646533333952701601768799639508076 ) ;
		m_rules[2][6] = std::make_pair( 6, 	.24261107190140773379964095790325635233 ) ;
	}	
}

