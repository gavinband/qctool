#ifndef TEXT_DIAGRAM_HPP
#define TEXT_DIAGRAM_HPP

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>

namespace textutils {
	struct Diagram
	{
	public:
		enum Border { Left = 0, Top = 1 , Right = 2, Bottom = 3 } ;
		enum BorderBehaviour { Clip, Clamp } ;
		enum PointType { Normal = '*', Dot = '.' } ;
		typedef std::pair< std::size_t, std::size_t > Dimensions ;
		typedef std::map< std::pair< double, double >, PointType > Points ;
		typedef std::map< std::pair< double, double >, PointType >::const_iterator PointsIterator ;

	public:
		Diagram()
		: m_border_behaviours(4, Clip)
			m_x_lower_limit( std::numeric_limits< double >::infinity()),
			m_y_lower_limit( std::numeric_limits< double >::infinity())
		{}
		
		virtual ~Diagram() {}
		
		void set_range( double x_lower_limit, double x_upper_limit, double y_lower_limit, double y_upper_limit ) {
			set_x_range( x_lower_limit, x_upper_limit ) ;
			set_y_range( y_lower_limit, y_upper_limit ) ;
		} ;

		void set_x_range( double x_lower_limit, double x_upper_limit ) {
			m_x_lower_limit = x_lower_limit ;
			m_x_upper_limit = x_upper_limit ;
		} ;

		void set_y_range( double y_lower_limit, double y_upper_limit ) {
			m_y_lower_limit = y_lower_limit ;
			m_y_upper_limit = y_upper_limit ;
		} ;

		void set_outer_dimensions( std::size_t width, std::size_t height ) {
			assert( width > 0 && height > 0 ) ;
			m_outer_dimensions = Dimensions( width, height ) ;
		}

		void set_border_behaviour( Border const border, BorderBehaviour const behaviour ) {
			m_border_behaviours[ border ] = behaviour ;
		}

		void clear() {
			m_points.resize(0) ;
		}

		void add_point( std::size_t const layer, double const x, double const y, PointType const type = Normal ) {
			if( !x_limits_have_been_set ) {
				m_x_lower_limit = std::min( m_x_lower_limit, x ) ;
				m_x_upper_limit = std::max( m_x_upper_limit, x ) ;
			}
			if( !y_limits_have_been_set ) {
				m_y_lower_limit = std::min( m_y_lower_limit, y ) ;
				m_y_upper_limit = std::max( m_y_upper_limit, y ) ;
			}
			m_points[ layer ][ std::make_pair( x, y ) ] = type ;
		}

		virtual std::string draw() const = 0 ;

		std::size_t add_layer() {
			m_points.push_back( Points() ) ;
			return m_points.size() - 1 ;
		}

	protected:
		
		std::pair< double, double > x_limits() const { return std::make_pair( m_x_lower_limit, m_x_upper_limit ) ; }
		std::pair< double, double > y_limits() const { return std::make_pair( m_y_lower_limit, m_y_upper_limit ) ; }
		Dimensions outer_dimensions() const { return m_outer_dimensions ; }

		std::vector< Points > const& points() const { return m_points ; }

		bool apply_clipping( double& x, double& y, Dimensions const& dimensions ) const {
			if( apply_1d_clipping( x, dimensions.first, Left, Right ) ) {
				return apply_1d_clipping( y, dimensions.second, Bottom, Top ) ;
			}
			return false ;
		}

	private:
		
		bool x_limits_have_been_set() const {
			return m_x_lower_limit != std::numeric_limits< double >::infinity() ;
		}

		bool y_limits_have_been_set() const {
			return m_y_lower_limit != std::numeric_limits< double >::infinity() ;
		}
		
		bool is_NaN( double const value ) const {
			return !(value == value );
		}
		
		bool apply_1d_clipping( double& value, double upper_limit, Border const lower_border, Border const upper_border ) const {
			if( is_NaN( value )) {
				return false ;
			}

			if( value < 0 ) {
				if( m_border_behaviours[ lower_border ] == Clamp ) {
					value = 0 ;
				}
				else {
					return false ;
				}
			}
			else if( value >= upper_limit ) {
				if( m_border_behaviours[ upper_border ] == Clamp ) {
					value = upper_limit - 1 ;
				}
				else {
					return false ;
				}
			}
			else {
				return true ;
			}
			return true ;
		}

	private:
	
		double m_x_lower_limit, m_x_upper_limit, m_y_lower_limit, m_y_upper_limit ;	
		Dimensions m_outer_dimensions ;
		std::vector< BorderBehaviour > m_border_behaviours ;
		std::vector< Points > m_points ;	
	} ;

	std::ostream& operator<<( std::ostream& out, Diagram const& diagram ) {
		return out << diagram.draw() ;
	}

	struct Canvas
	{
		virtual ~Canvas() {}
		virtual char& operator()( std::size_t i, std::size_t j ) = 0 ;
		virtual char const& operator()( std::size_t i, std::size_t j ) const = 0 ;

		virtual std:size_t width() const = 0 ;
		virtual std:size_t height() const = 0 ;
		
		std::string draw() const {
			std::string result ;
			for( std::size_t j = height(); j > 0 ; --j ) {
				for( std::size_t i = 0; i < width() ; ++i ) {
						result += canvas( i, j - 1 ) ;
				}
				result += "\n" ;
			}
			return result ;
		}

		
	} ;

	struct CanvasView: public Canvas
	{
	public:
		CanvasView( Canvas& underlying_canvas, std::size_t x_origin, std::size_t y_origin, std::size_t width, std::size_t height )
			: m_underlying_canvas( canvas ),
			m_x_origin( x_origin ),
			m_y_origin( y_origin ),
			m_width( width ),
			m_height( height )
		{
			assert( m_x_origin + m_width <= m_underlying_canvas.width() ) ;
			assert( m_y_origin + m_height <= m_underlying_canvas.height() ) ;
		}
			
		char& operator()( std::size_t i, std::size_t j ) {
			assert( i < m_width && j < m_height ) ;
			return m_underlying_canvas( i + m_x_origin, j + m_y_origin ) ;
		}
			
		char const& operator()( std::size_t i, std::size_t j ) const {
			assert( i < m_width && j < m_height ) ;
			return m_underlying_canvas( i + m_x_origin, j + m_y_origin ) ;
		}

		std:size_t width() const { return m_width ; }
		std:size_t height() const { return m_height ; }

	private:
			
		Canvas& m_underlying_canvas ;
		std::size_t
			m_x_origin,
			m_y_origin,
			m_width,
			m_height ;
	} ;

	struct VectorCanvas: public Canvas
	{
		VectorCanvas( std::size_t width, std::size_t height )
			: m_height( height )
		{
			m_data.resize( width ) ;
			for( std::size_t i = 0; i < width; ++i ) {
				m_data[i].resize( height ) ;
			}
		}
		
		char& operator()( std::size_t i, std::size_t j ) { return m_data[i][j] ; }
		char const& operator()( std::size_t i, std::size_t j ) const { return m_data[i][j] ; }

		std::size_t width() const { return m_data.size() ; }
		std::size_t height() const { return m_height ; }

	private:
		std::size_t m_height ;
		std::vector< std::vector< char > > m_data ;
	} ;


	/*
	* A Diagram with axes, looking something like this:
	*
	* 10 |
	*    |              +
	*    |     +++   + 
	*    |                      +
	* 0  +------------------------------------
	*    0                                 100
	*/
	struct SimpleDiagramWithAxes: public Diagram
	{
	public:

		SimpleDiagramWithAxes( Canvas* external_canvas )
			: m_external_canvas( external_canvas )
		{
			assert( m_external_canvas != 0 ) ;
			set_outer_dimensions( external_canvas->width(), external_canvas->height() ) ;
		}

		SimpleDiagramWithAxes( std::size_t width, std::size_t height )
			: m_external_canvas( 0 )
		{
			m_internal_canvas.resize( width, height ) ;
			set_outer_dimensions( width, height ) ;
		}

		Canvas const& draw() const {
			Dimensions canvas_dimensions = outer_dimensions() ;
			Canvas& canvas = get_canvas() ;
			Dimensions axes_dimensions = get_space_needed_for_axes() ;
			draw_axes( axes_dimensions, canvas_dimensions, canvas ) ;
			draw_points( axes_dimensions, canvas_dimensions, canvas ) ;
			return canvas ;
		}
		
	private:
		Canvas& get_canvas() const {
			if( m_external_canvas != 0 ) {
				return *m_external_canvas ;
			}
			else {
				return m_internal_canvas ;
			}
		}
		
		Dimensions get_space_needed_for_axes() const {
			std::string y_lower_limit = to_string( y_limits().first ) ;
			std::string y_upper_limit = to_string( y_limits().second ) ;
			return Dimensions(
				std::max( y_lower_limit.size(), y_upper_limit.size() ) + 2,
				2
			) ;
		}
		
		std::string to_string( double const value ) const {
			std::ostringstream ostr ;
			ostr << std::fixed << value ;
			return ostr.str() ;
		}
		
		void draw_axes( Dimensions const& axes_dimensions, Dimensions const& canvas_dimensions, Canvas& canvas ) const {
			std::cout << axes_dimensions.first << ":" << axes_dimensions.second << "\n" ;
			canvas( axes_dimensions.first - 1, axes_dimensions.second - 1 ) = '+' ;
			for( std::size_t i = axes_dimensions.first; i < canvas_dimensions.first; ++i ) {
				canvas( i, axes_dimensions.second - 1 ) = '-' ;
			}
			for( std::size_t i = axes_dimensions.second; i < canvas_dimensions.second; ++i ) {
				canvas( axes_dimensions.first - 1, i ) = '|' ;
			}
			draw_x_limits_if_there_is_space( axes_dimensions, canvas_dimensions, canvas ) ;
			draw_y_limits_if_there_is_space( axes_dimensions, canvas_dimensions, canvas ) ;
		}
		
		void draw_x_limits_if_there_is_space( Dimensions const& axes_dimensions, Dimensions const& canvas_dimensions, Canvas& canvas ) const {
			std::string x_lower_limit = to_string( x_limits().first ) ;
			std::string x_upper_limit = to_string( x_limits().second ) ;
			
			if( canvas_dimensions.first > ( x_lower_limit.size() + axes_dimensions.first - 1 )) {
				for( std::size_t i = 0; i < x_lower_limit.size(); ++i ) {
					canvas( axes_dimensions.first - 1 + i, 0 ) = x_lower_limit[i] ;
				}
			}
			
			if( canvas_dimensions.first > ( x_lower_limit.size() + x_upper_limit.size() + axes_dimensions.first - 1 + 1 )) {
				for( std::size_t i = 0; i < x_upper_limit.size(); ++i ) {
					canvas( canvas_dimensions.first - x_upper_limit.size() + i, 0 ) = x_upper_limit[i] ;
				}
			}
		}

		void draw_y_limits_if_there_is_space( Dimensions const& axes_dimensions, Dimensions const& canvas_dimensions, Canvas& canvas ) const {
			std::string y_lower_limit = to_string( y_limits().first ) ;
			std::string y_upper_limit = to_string( y_limits().second ) ;
			
			if( canvas_dimensions.second >= axes_dimensions.second ) {
				for( std::size_t i = 0; i < y_lower_limit.size(); ++i ) {
					canvas( i, axes_dimensions.second - 1 ) = y_lower_limit[i] ;
				}
			}
			
			if( canvas_dimensions.second > axes_dimensions.second + 1 ) {
				for( std::size_t i = 0; i < y_upper_limit.size(); ++i ) {
					canvas( i, canvas_dimensions.second - 1 ) = y_upper_limit[i] ;
				}
			}
		}


		void draw_points( Dimensions const& axes_dimensions, Dimensions const& canvas_dimensions, Canvas& canvas) const {
			for( std::vector<Points>::const_iterator i = points().begin() ; i != points().end(); ++i ) {
				draw_points( *i, axes_dimensions, canvas_dimensions, canvas ) ;
			}
		}
		
		void draw_points( Points const& points, Dimensions const& axes_dimensions, Dimensions const& canvas_dimensions, Canvas& canvas) const {
			Dimensions inner_dimensions (
				canvas_dimensions.first - axes_dimensions.first,
				canvas_dimensions.second - axes_dimensions.second
			) ;
			
			double x_scale = inner_dimensions.first / ( x_limits().second - x_limits().first ) ;
			double y_scale = inner_dimensions.second / ( y_limits().second - y_limits().first ) ;

			for( Points::const_iterator i = points.begin() ; i != points.end(); ++i ) {
				double
					x = (i->first.first - x_limits().first) * x_scale,
					y = (i->first.second - y_limits().first) * y_scale ;
					
				if( apply_clipping( x, y, inner_dimensions )) {
					canvas( axes_dimensions.first + x, axes_dimensions.second + y ) = i->second ;
				}
			}
		}
		
		Canvas* m_external_canvas ;
		VectorCanvas m_internal_canvas ;
	} ;
}


#endif
