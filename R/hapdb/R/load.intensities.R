load.intensities <-
function( hapdb, chromosome = NULL, rsids = NULL, range = NULL, positions = NULL, samples = NULL, analysis = NULL, verbose = FALSE ) {
	require( RSQLite )
	require( Rcompression )
	sql = paste(
		"SELECT analysis_id, E.name AS analysis, V.id AS variant_id, chromosome, position, rsid, alleleA, alleleB, H.N AS N, H.data",
		"FROM Entity E",
		"INNER JOIN Variant V",
		"CROSS JOIN Intensity H ON H.analysis_id == E.id AND H.variant_id == V.id",
		sep = " "
	)
	
	if( is.null( analysis ) ) {
		analysis = as.character( unique( hapdb$samples$analysis ) )
		if( length( analysis ) > 1 ) {
			stop( "This hapdb file has more than one analysis, please specify one." ) ;
		} else if( length( analysis ) == 0 ) {
			stop( "This hapdb file seems to have no analyses!" ) ;
		}
	}
	
	analysis_id = hapdb$samples$analysis_id[ match( analysis, hapdb$samples$analysis ) ]
	sql = paste(
		sql,
		sprintf( "WHERE analysis_id == %d", analysis_id ),
		sep = " "
	) ;
   
	if(( !is.null( chromosome ) || !is.null( range ) ) && !is.null( positions ) ) {
        stop( "You can't specify chromosome and range as well as a list of positions." )
    }
	if( !is.null( chromosome ) ) {
		sql = paste(
			sql,
			sprintf( "AND chromosome == '%s'", as.character( chromosome ) ),
			sep = " "
		)
	}
	if( !is.null( range ) ) {
		sql = paste(
			sql,
			sprintf( "AND position BETWEEN %d AND %d", as.integer( range[1] ), as.integer( range[2] ) ),
			sep = " "
		)
	}
	if( !is.null( rsids ) ) {
		sql = paste(
			sql,
			"AND rsid IN (",
			paste( sprintf( "'%s'", rsids ), collapse = ',' ),
			')',
			sep = " "
		)
	}
    if( !is.null( positions )) {
        sql = paste(
            sql,
            "AND ( ",
            sep = " " 
        )   
        for( i in 1:nrow( positions ) ) { 
            if( i > 1 ) { 
                sql = paste( sql, "OR", sep = " " )
            }   
            sql = paste(
                sql,
                sprintf( "( chromosome = '%s' AND position == '%d' )", positions[i,1], positions[i,2] ),
                sep = " " 
            )   
        }   
        sql = paste(
            sql,
            ")",
            sep = " " 
        )   
    }   
	if( verbose ) {
		cat( "load.intensities(): running query :\"", sql, "\"...\n", sep = "" ) ;
		print( dbGetQuery( hapdb$db, sprintf( "EXPLAIN QUERY PLAN %s", sql ) ) )
	}
	D = dbGetQuery( hapdb$db, sql )

	if( nrow(D) > 0 ) {
    	if( verbose ) {
    		cat( "load.intensities(): loaded, uncompressing...\n" ) ;
        }
		result = list(
			variant = D[,-which( colnames(D) == "data")],
			data = matrix( NA, nrow = nrow(D), ncol = 2 * D$N[1] )
		)

		for( i in 1:nrow(D) ) {
			compressed_data = unlist( D$data[i] )
			uncompressed_data = uncompress( compressed_data, asText = FALSE )
			result$data[i,] = parse_intensities( uncompressed_data, D$N[i] )
		}
	} else {
		result = list(
			variant = D[,-which( colnames(D) == "data")],
			data = matrix( NA, nrow = 0, ncol = 2 * length( which( hapdb$samples$analysis == analysis ) ) )
		) ;
	}
	result$samples = hapdb$samples[ which( hapdb$samples$analysis == analysis ), ]
	if( !is.null( samples ) ) {
		if( mode( samples ) == "character" ) {
			samples.choice = sort( which( hapdb$samples$analysis == analysis & hapdb$samples$identifier %in% samples ) ) ;
		} else {
			samples.choice = samples
		}
		genotype.choice = sort( union( (samples.choice * 2) - 1, samples.choice * 2 ) )
	
		result$data = result$data[, genotype.choice ]
		result$samples = result$samples[ samples.choice, ]
		if( nrow( result$variant ) > 0 ) {
			result$variant$N = length( samples )
		}
	}

	colnames( result$data ) = rep( NA, 2 * nrow( result$sample ) )
	colnames( result$data )[ seq( from = 1, by = 2, length = nrow( result$samples ) ) ] = paste( result$samples[, 'identifier' ], "X", sep = ":" )
	colnames( result$data )[ seq( from = 2, by = 2, length = nrow( result$samples ) ) ] = paste( result$samples[, 'identifier' ], "Y", sep = ":" )

	if( verbose ) {
		cat( "load.intensities(): done.\n" ) ;
	}
	return( result )
}

