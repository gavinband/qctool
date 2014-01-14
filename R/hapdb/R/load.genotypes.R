load.genotypes <-
function( hapdb, chromosome = NULL, positions = NULL, rsid = NULL, range = NULL, samples = NULL, analysis = NULL, verbose = FALSE, dosage.threshhold = NULL ) {
	require( RSQLite )
	require( Rcompression )
	sql = paste(
		"SELECT analysis_id, E.name AS analysis, V.id AS variant_id, chromosome, position, rsid, alleleA, alleleB, H.N AS N, H.data",
		"FROM Entity E",
		"INNER JOIN Variant V",
		"CROSS JOIN Genotype H ON H.analysis_id == E.id AND H.variant_id == V.id",
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
	if( !is.null( rsid ) ) {
		sql = paste(
			sql,
			"AND rsid IN (",
			paste( sprintf( "'%s'", rsid ), collapse = ',' ),
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
		cat( "load.genotypes(): running query :\"", sql, "\"...\n", sep = "" ) ;
		print( dbGetQuery( hapdb$db, sprintf( "EXPLAIN QUERY PLAN %s", sql ) ) )
	}
	D = dbGetQuery( hapdb$db, sql )

	if( nrow(D) > 0 ) {
    	if( verbose ) {
    		cat( "load.genotypes(): loaded, uncompressing...\n" ) ;
        }
		result = list(
			variant = D[,-which( colnames(D) == "data")],
			data = matrix( NA, nrow = nrow(D), ncol = 3 * D$N[1] )
		)

		for( i in 1:nrow(D) ) {
			compressed_data = unlist( D$data[i] )
			uncompressed_data = uncompress( compressed_data, asText = FALSE )
			result$data[i,] = parse_genotypes( uncompressed_data, D$N[i] )
		}
	} else {
		result = list(
			variant = D[,-which( colnames(D) == "data")],
			data = matrix( NA, nrow = 0, ncol = 3 * length( which( hapdb$samples$analysis == analysis ) ) )
		) ;
	}
	result$samples = hapdb$samples[ which( hapdb$samples$analysis == analysis ), ]
	if( !is.null( samples ) ) {
		if( mode( samples ) == "character" ) {
			samples.choice = sort( which( hapdb$samples$analysis == analysis & hapdb$samples$identifier %in% samples ) ) ;
		} else {
			samples.choice = samples
		}
		genotype.choice = sort( union( (samples.choice * 3) - 2, (samples.choice * 3) - 1, samples.choice * 3 ) )
	
		result$data = result$data[, genotype.choice ]
		result$samples = result$samples[ samples.choice, ]
		if( nrow( result$variant ) > 0 ) {
			result$variant$N = length( samples )
		}
	}
	colnames( result$data ) = rep( NA, 3 * nrow( result$sample ) )
	colnames( result$data )[ seq( from = 1, by = 3, length = nrow( result$samples ) ) ] = paste( result$samples[, 'identifier' ], "AA", sep = ":" )
	colnames( result$data )[ seq( from = 2, by = 3, length = nrow( result$samples ) ) ] = paste( result$samples[, 'identifier' ], "AB", sep = ":" )
	colnames( result$data )[ seq( from = 3, by = 3, length = nrow( result$samples ) ) ] = paste( result$samples[, 'identifier' ], "BB", sep = ":" )

	for( i in 1:nrow( result$samples ) ) {
		# handle the case of three zeroes.
		# this shouldn't occur with the latest qctool, but did in earlier versions (as it does in bgen etc.)
		w = which( rowSums( result$data[, ((3*i)-2):(3*i)] ) == 0 )	
		if( length(w) > 0 ) {
			result$data[ w, ((3*i)-2):(3*i)] = NA
		}
	}
	
	if( !is.null( dosage.threshhold )) {
	    N = nrow( result$samples )
	    dosage = matrix( NA, nrow = nrow( result$data ), ncol = ncol( result$data ) / 3 ) ;
        A = result$data ;
        A[ which( A >= dosage.threshhold )] = 1 ;
        A[ which( A < dosage.threshhold )] = 0 ;
	    dosage = A[ , seq( from = 2, by = 3, length = N ) ] + 2 * A[ , seq( from = 3, by = 3, length = N ) ] ;
	    result$dosage = dosage ;
	    colnames( result$dosage ) = result$samples$identifier
	}
	
	if( verbose ) {
		cat( "load.genotypes(): done.\n" ) ;
	}
	return( result )
}

