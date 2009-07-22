#ifndef CASE_CONTROL_STATISTIC_HPP
#define CASE_CONTROL_STATISTIC_HPP

// Base class for individual statistics
struct CaseControlStatistic
{
	public:
		CaseControlStatistic() {} ;
		virtual ~CaseControlStatistic() {}

		typedef std::map< double, GenotypeAmounts > case_status_statistic_map_t ;

		template< typename T >
		T get_value( case_status_statistic_map_t const& ) const ;

	protected:
		virtual double calculate_value( case_status_statistic_map_t const& ) const = 0;
} ;


template<>
double CaseControlStatistic::get_value< double >( case_status_statistic_map_t const& case_status_statistic_map ) const
{
	return calculate_value( case_status_statistic_map ) ;
}

#endif
