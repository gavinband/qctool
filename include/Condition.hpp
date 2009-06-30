
#ifndef __GTOOL_CONDITION_HPP__
#define __GTOOL_CONDITION_HPP__

#include <set>
#include <vector>
#include "GToolException.hpp"

struct ConditionException: public GToolException
{
	ConditionException( std::string const& msg ) ;
} ;

// Condition is a base class for a polymorphic hierarchy of tests to be performed on some data.
// The data is broken up into mandatory, and optional data.  This is intended as a relatively flexible mechanism
// to allow tests to share calculated data (but not require it if not needed).
template< typename MandatoryData, typename OptionalData >
class Condition
{
	public:
		Condition() {} ;
		virtual ~Condition() {} ;

		virtual bool check_if_satisfied( MandatoryData const&, OptionalData const* = 0 ) const = 0;

		virtual void format_to_stream( std::ostream& ) const = 0 ;

	private:
	
		// prevent copying and assignment.
		Condition( Condition const& ) ;
		Condition& operator=( Condition const& ) ;
} ;

template< typename MandatoryData, typename OptionalData >
std::ostream& operator<< ( std::ostream& oStream, Condition< MandatoryData, OptionalData > const& condition ) {
	condition.format_to_stream( oStream ) ;
	return oStream ;
} 


template< typename MandatoryData, typename OptionalData >
struct InvertCondition: public Condition< MandatoryData, OptionalData >
{
	typedef Condition< MandatoryData, OptionalData > base_t ;

    InvertCondition( std::auto_ptr< base_t > condition_to_invert )
		: m_condition_to_invert( condition_to_invert )
	{
		assert( m_condition_to_invert.get() != 0 ) ;
	}

	bool check_if_satisfied( MandatoryData const& mandatory_data , OptionalData const * optional_data = 0 ) const {
		return !( m_condition_to_invert->check_if_satisfied( mandatory_data, optional_data )) ;
	}

	void format_to_stream( std::ostream& oStream ) const {
		oStream
			<< "NOT "
			<< (*m_condition_to_invert) ;
	}

    private:
        std::auto_ptr< base_t > m_condition_to_invert ;
} ;


// Base-class for compound conditions such as implementations of logical AND, or OR.
template< typename MandatoryData, typename OptionalData >
struct CompoundCondition: public Condition< MandatoryData, OptionalData >
{
	public:
		typedef Condition< MandatoryData, OptionalData > base_t ;
		typedef typename std::vector< base_t* >::const_iterator subcondition_iterator_t ;

	public:
		~CompoundCondition() {
			typename std::vector< base_t* >::iterator
				subcondition_i( m_subconditions.begin() ),
				subcondition_end( m_subconditions.end() ) ;
				
			for( ; subcondition_i != subcondition_end; ++subcondition_i ) {
				delete (*subcondition_i) ;
			}
		}

		void add_subcondition( std::auto_ptr< Condition< MandatoryData, OptionalData > > subcondition ) {
			assert( subcondition.get() != 0 ) ;
			m_subconditions.push_back( subcondition.release() ) ;
		}

		void format_to_stream( std::ostream& oStream, std::string const& join ) const {
			subcondition_iterator_t
				first_subcondition_i( begin_subconditions() ),
				subcondition_i( begin_subconditions() ),
				subcondition_end( end_subconditions() ) ;

			for( ; subcondition_i != subcondition_end; ++subcondition_i ) {
				if( subcondition_i != first_subcondition_i ) {
					oStream << " " << join << " " ;
				}
				oStream << (**subcondition_i) ;
			}
		}


		subcondition_iterator_t begin_subconditions() const { return m_subconditions.begin() ; }
		subcondition_iterator_t end_subconditions() const { return m_subconditions.end() ; }

	private:
	
		std::vector< base_t* > m_subconditions ;
} ;


template< typename MandatoryData, typename OptionalData >
struct AndCondition: public CompoundCondition< MandatoryData, OptionalData >
{
	typedef CompoundCondition< MandatoryData, OptionalData > base_t ;
	
	bool check_if_satisfied( MandatoryData const& mandatory_data, OptionalData const * optional_data = 0 ) const {
		for( typename base_t::subcondition_iterator_t subcondition_i = base_t::begin_subconditions() ; subcondition_i != base_t::end_subconditions(); ++subcondition_i ) {
			if( !(*subcondition_i)->check_if_satisfied( mandatory_data, optional_data ) ) {
				return false ;
			}
		}
		return true ;
	}
	
	void format_to_stream( std::ostream& oStream ) const {
		base_t::format_to_stream( oStream, "AND" ) ;
	}
} ;


template< typename MandatoryData, typename OptionalData >
struct OrCondition: public CompoundCondition< MandatoryData, OptionalData >
{
	typedef CompoundCondition< MandatoryData, OptionalData > base_t ;

	bool check_if_satisfied( MandatoryData const& mandatory_data, OptionalData const * optional_data = 0 ) const {
		bool result = false ;
		for( typename base_t::subcondition_iterator_t subcondition_i = base_t::begin_subconditions() ; subcondition_i != base_t::end_subconditions(); ++subcondition_i ) {
			result = result || (*subcondition_i)->check_if_satisfied( mandatory_data, optional_data ) ;
		}
		return result ;
	}
	
	void format_to_stream( std::ostream& oStream ) const {
		base_t::format_to_stream( oStream, "OR" ) ;
	}
} ;


#endif

