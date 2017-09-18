$(document).ready(
	function() {
		$( '.nav_button' ).hover(
			function( elt ) {
				$(elt.currentTarget).children( 'ul' ).addClass( "open" ) ;
			},
			function( elt ) {
				$(elt.currentTarget).children( 'ul' ).removeClass( "open" ) ;
			}
		)
	}
) ;
