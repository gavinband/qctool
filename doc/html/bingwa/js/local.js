var globals = {
	current_page: "none"
} ;

var setPage = function( page ) {
	$( ".nav_button" ).children( "span" ).each( function( i ) { $(this).removeClass( "nav_button_selected" ) ; } ) ;
	$( '[name=' + page + '].nav_button' ).children( 'span' ).addClass( "nav_button_selected" ) ;
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
	}
) ;
