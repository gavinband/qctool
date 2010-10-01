var globals = {
	current_page: "none"
} ;

var setPage = function( page ) {
	$( ".nav_button" ).each( function( i ) { $(this).css( "border-bottom", "none" ) ; } ) ;
	$( '[name=' + page + '].nav_button' ).css( "border-bottom", "2px solid black" ) ;
	$( '#' + page ).fadeIn( 30 ) ;
}

var changePage = function( page ) {
	if( page != globals.current_page ) {
		if( globals.current_page != "none" ) {
			$( '#' + globals.current_page ).fadeOut( 30,
				function() {
					setPage( page ) ;
				} ) ;
		} else {
			setPage( page ) ;
		}
		globals.current_page = page ;
		location.hash = page ;
	}
}

var setFirstPage = function( default_page ) {
	if( location.hash == "" ) {
		location.hash = "#" + default_page ;
	} else {
		$( window ).trigger( 'hashchange' ) ;
	}
}

$(document).ready(
	function() {
		$( ".nav_button" ).each( function( i ) { $(this).attr( "inactive_src", "style/images/" + $(this).attr( "name" ) + ".png" ) } )
		$( ".nav_button" ).each( function( i ) { $(this).attr( "active_src", "style/images/" + $(this).attr( "name" ) + "_red.png" ) } )
		
		$( ".nav_button" ).hover(
			function( eventObject ) {
				name = $(this).attr( "name" ) ;
				$(this).attr( "src", $(this).attr( "active_src" ) ) ;
			},
			function( eventObject ) {
				name = $(this).attr( "name" ) ;
				$(this).attr( "src", $(this).attr( "inactive_src" )) ;
			}
		)
		
		$( ".nav_button" ).click(
			function( eventObject ) {
				page = $( this ).attr( 'name' ) ;
				location.hash = page ;
			}
		) ;
		
		$(window).bind( 'hashchange', function() {
			// Get the page.  Omit the leading '#'.
			var page = location.hash.substring( 1 ) ;
			changePage( page ) ;
		} ) ;
		
		// Setup up the resizeable terminal-display things
		$(  ".terminal_display" ).attr( "display_size", "small" ) ;
		$( ".terminal_display" ).hover(
			function( eventObject ) {
				$( this ).css( "background-color", "#10A145" ) ;
			},
			function( eventObject ) {
				$( this ).css( "background-color", "#22AD57" ) ;
			}
		) ;
		$( ".terminal_display" ).click(
			function( eventObject ) {
				switch_value = $(this).attr( "display_size" ) ;
				if( switch_value == "small" ) {
					$( this ).animate(
						{ height: "500px" },
						duration = 100,
						callback = function() {
							$(this).attr( "display_size", "large" ) ;
						}
					) ;
				}
				else {
					$( this ).animate(
						{ height: "110px" },
						duration = 100,
						callback = function() {
							$(this).attr( "display_size", "small" ) ;
						}
					) ;
				}
			}
		) ;
	}
) ;
