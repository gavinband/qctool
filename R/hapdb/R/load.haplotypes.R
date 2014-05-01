load.haplotypes <-
function(
	hapdb,
	chromosome = NULL, range = NULL, rsids = NULL, positions = NULL, variant_ids = NULL,
	samples = NULL, analysis = NULL,
	verbose = FALSE,
	compute.dosage = FALSE,
	method = "cpp"
) {
	require( RSQLite )
	require( Rcompression )
	sql = paste(
		"SELECT analysis_id, E.name AS analysis, V.id AS variant_id, chromosome, position, rsid, alleleA, alleleB, H.N AS N, H.data",
		"FROM Entity E",
		"INNER JOIN Variant V",
		"CROSS JOIN Haplotype H ON H.analysis_id == E.id AND H.variant_id == V.id",
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
	
	dbGetQuery( hapdb$db, "CREATE TEMPORARY TABLE tmpHapdbLoadGenotypes ( variant_id INT NOT NULL )" ) ;
	if( !is.null( chromosome ) && !is.null( range )) {
		dbGetQuery( hapdb$db, sprintf( "INSERT INTO tmpHapdbLoadGenotypes SELECT id FROM Variant WHERE chromosome == '%s' AND position BETWEEN %d AND %d", chromosome, range[1], range[2] ) ) ;
    }
	if( !is.null( rsids ) ) {
		dbGetPreparedQuery( hapdb$db, "INSERT INTO tmpHapdbLoadGenotypes SELECT id FROM Variant WHERE rsid == ?", rsids ) ;
	}
    if( !is.null( positions )) {
		dbGetPreparedQuery( hapdb$db, "INSERT INTO tmpHapdbLoadGenotypes SELECT id FROM Variant WHERE chromosome == ? AND position == ?", positions ) ;
    }
    if( !is.null( variant_ids )) {
		dbGetPreparedQuery( hapdb$db, "INSERT INTO tmpHapdbLoadGenotypes VALUES(?)", data.frame( variant_id = variant_ids ) ) ;
    }
    sql = paste(
        sql,
        "AND V.id IN ( SELECT variant_id FROM tmpHapdbLoadGenotypes )",
        sep = " "
    )

	if( verbose ) {
		cat( "load.haplotypes(): running query :\"", sql, "\"...\n", sep = "" ) ;
		print( dbGetQuery( hapdb$db, sprintf( "EXPLAIN QUERY PLAN %s", sql ) ) )
	}

	#############
	# Get the data
	D = dbGetQuery( hapdb$db, sql )
	dbGetQuery( hapdb$db, "DROP TABLE tmpHapdbLoadGenotypes" ) ;

	#############
	# Now post-process
	if( nrow(D) > 0 ) {
		result = list(
			variant = D[,-which( colnames(D) == "data"), drop = FALSE]
		)

		if( method == "cpp" ) {
			if( verbose ) {
				cat( "load.haplotypes(): uncompressing using C++...\n" ) ;
			}
		    result$data = rcpp_uncompress_bitpack_haplotypes( D$data, D$N[1] ) ;
    	} else {
			if( verbose ) {
			    cat( "load.haplotypes(): uncompressing using R...\n" ) ;
			}
			result$data = matrix( NA, nrow = nrow(D), ncol = 2 * D$N[1] )
    		for( i in 1:nrow(D) ) {
    			compressed_data = unlist( D$data[i] )
    			uncompressed_data = uncompress( compressed_data, asText = FALSE )
    			result$data[i,] = parse_haplotypes( uncompressed_data, D$N[i] )
    			if( i %% 100 == 0 ) {
    				gc()
    				if( verbose ) {
    					cat( "(done ", i, " of ", nrow(D), ")...\n", sep = "" ) ;
    				}
    			}
    		}
    	}
	} else {
		result = list(
			variant = D[,-which( colnames(D) == "data"), drop = FALSE],
			data = matrix( NA, nrow = 0, ncol = 2 * length( which( hapdb$samples$analysis == analysis ) ) )
		) ;
	}
	result$samples = hapdb$samples[ which( hapdb$samples$analysis == analysis ), ]
	N = nrow( result$samples )
	if( !is.null( samples ) ) {
		if( mode( samples ) == "character" ) {
			samples.choice = sort( which( hapdb$samples$analysis == analysis & hapdb$samples$identifier %in% samples ) ) ;
		} else {
			samples.choice = samples
		}
		hap.choice = sort( union( (samples.choice * 2) - 1, samples.choice * 2 ) )
	
		result$data = result$data[, hap.choice, drop = FALSE ]
		result$samples = result$samples[ samples.choice,, drop = FALSE ]
		if( nrow( result$variant ) > 0 ) {
			result$variant$N = length( samples )
		}
	}

    N = nrow( result$samples )
	colnames( result$data ) = rep( NA, 2 * nrow( result$sample ) )
	colnames( result$data )[ seq( from = 1, by = 2, length = N ) ] = paste( result$samples[, 'identifier' ], "0", sep = ":" )
	colnames( result$data )[ seq( from = 2, by = 2, length = N ) ] = paste( result$samples[, 'identifier' ], "1", sep = ":" )

	if( compute.dosage ) {
		A = result$data[, seq( from = 1, by = 2, length = N ), drop = FALSE ]
		B = result$data[, seq( from = 2, by = 2, length = N ), drop = FALSE ]
		B[ which( is.na( B ) ) ] = 0
		dosage = A + B
        colnames( dosage ) = result$samples$identifier
        rownames( dosage ) = result$variant$rsid
	    result$dosage = dosage
	}

	if( verbose ) {
		cat( "load.haplotypes(): done.\n" ) ;
	}
	return( result )
}

