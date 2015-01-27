bgen.read_snp <-
function( header ) {
 data = bgen.read_snp_id_data( header )
 if( is.null( data )) {
     return( NULL ) ;
 }
 else {
     return( bgen.read_snp_probability_data( header, data )) ;
 }
}
