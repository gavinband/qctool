
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef __GTOOL_CONDITION_HPP__
#define __GTOOL_CONDITION_HPP__

#include <set>
#include <vector>
#include <memory>
#include "GToolException.hpp"

struct ConditionException: public GToolException
{
	ConditionException( std::string const& msg ) ;
} ;

// Condition is a base class for a polymorphic hierarchy of tests to be performed on some data.
// The data is broken up into mandatory, and optional data.  This is intended as a relatively flexible mechanism
// to allow tests to share calculated data (but not require it if not needed).
template< typename MandatoryData >
class Condition
{
	public:
		Condition() {} ;
		virtual ~Condition() {} ;

		virtual bool check_if_satisfied( MandatoryData const& ) const = 0;

		virtual void format_to_stream( std::ostream& ) const = 0 ;

	private:
	
		// prevent copying and assignment.
		Condition( Condition const& ) ;
		Condition& operator=( Condition const& ) ;
} ;

template< typename MandatoryData >
std::ostream& operator<< ( std::ostream& oStream, Condition< MandatoryData > const& condition ) {
	condition.format_to_stream( oStream ) ;
	return oStream ;
} 


template< typename MandatoryData >
struct InvertCondition: public Condition< MandatoryData >
{
	typedef Condition< MandatoryData > base_t ;

    InvertCondition( std::auto_ptr< base_t > condition_to_invert )
		: m_condition_to_invert( condition_to_invert )
	{
		assert( m_condition_to_invert.get() != 0 ) ;
	}

	bool check_if_satisfied( MandatoryData const& mandatory_data ) const {
		return !( m_condition_to_invert->check_if_satisfied( mandatory_data )) ;
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
template< typename MandatoryData >
struct CompoundCondition: public Condition< MandatoryData >
{
	public:
		typedef Condition< MandatoryData > base_t ;
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

		void add_subcondition( std::auto_ptr< Condition< MandatoryData > > subcondition ) {
			assert( subcondition.get() != 0 ) ;
			m_subconditions.push_back( subcondition.release() ) ;
		}

		void format_to_stream( std::ostream& oStream, std::string const& join ) const {
			if( m_subconditions.empty() ) {
				oStream << "(none)" ;
			}
			else {
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
		}

		std::size_t number_of_subconditions() const { return m_subconditions.size() ; }
		base_t& subcondition( std::size_t i ) const { return *m_subconditions[i] ; }
		subcondition_iterator_t begin_subconditions() const { return m_subconditions.begin() ; }
		subcondition_iterator_t end_subconditions() const { return m_subconditions.end() ; }

	private:
	
		std::vector< base_t* > m_subconditions ;
} ;


template< typename MandatoryData >
struct AndCondition: public CompoundCondition< MandatoryData >
{
	typedef CompoundCondition< MandatoryData > base_t ;
	
	bool check_if_satisfied( MandatoryData const& mandatory_data ) const {
		for( typename base_t::subcondition_iterator_t subcondition_i = base_t::begin_subconditions() ; subcondition_i != base_t::end_subconditions(); ++subcondition_i ) {
			if( !(*subcondition_i)->check_if_satisfied( mandatory_data ) ) {
				return false ;
			}
		}
		return true ;
	}
	
	void format_to_stream( std::ostream& oStream ) const {
		base_t::format_to_stream( oStream, "AND" ) ;
	}
} ;


template< typename MandatoryData >
struct OrCondition: public CompoundCondition< MandatoryData >
{
	typedef CompoundCondition< MandatoryData > base_t ;

	bool check_if_satisfied( MandatoryData const& mandatory_data ) const {
		bool result = false ;
		for( typename base_t::subcondition_iterator_t subcondition_i = base_t::begin_subconditions() ; subcondition_i != base_t::end_subconditions(); ++subcondition_i ) {
			result = result || (*subcondition_i)->check_if_satisfied( mandatory_data ) ;
		}
		return result ;
	}
	
	void format_to_stream( std::ostream& oStream ) const {
		base_t::format_to_stream( oStream, "OR" ) ;
	}
} ;


#endif

