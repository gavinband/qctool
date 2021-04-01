$(document).ready(
	function() {
		$( '.nav_button' ).hover(
			function( elt ) {
				$(elt.currentTarget).addClass( "open" ) ;
			},
			function( elt ) {
				$(elt.currentTarget).removeClass( "open" ) ;
			}
		)
	}
) ;
