library( GWAS )
cohort1 = load.sample.file( "~/Projects/Software/marchini/snptest/branches/v2.3.0/example/cohort1.sample")
cohort2 = load.sample.file( "~/Projects/Software/marchini/snptest/branches/v2.3.0/example/cohort1.sample")
cohort = rbind( cohort1, cohort2 )
cohort$sex = sample( size = nrow( cohort ), x = c( "male", "female" ), replace = T, prob = c(0.5, 0.5) )

compute.snp_stats <- function(
	snps,
	probs, # a matrix of probs, SNPs on rows, samples having three columns each.
	sex
) {
	result = snps ;
	result$a_allele_freq = NA
	result$b_allele_freq = NA
	result$missing_proportion = NA

	stopifnot( nrow( snps ) == nrow( probs )) ;
	for( i in 1:nrow( snps )) {
		if( (i-1) %% 100 == 0 ) {
			cat( i, "of", nrow(snps), "...\n" )
		}
		G = matrix( data = probs[i,], byrow = T, ncol = 3 ) ;

		if( snps$chromosome == "0X" ) {
			G.male = G
			G.male[ which( sex != "male" ), ] = 0
			G.female = G
			G.female[ which( sex != "female" ), ] = 0
		
			b_allele_count = sum( G.male[,2] ) + sum( G.female[,2] ) + 2 * sum( G.female[,3] )
			b_allele_freq = b_allele_count / ( sum( G.male ) + 2 * sum( G.female ))
			a_allele_freq = 1 - b_allele_freq ;
		} else {
			b_allele_count = sum( G[,2] ) + 2 * sum( G[,3] )
			b_allele_freq = b_allele_count / ( 2 * sum( G ) )
			a_allele_freq = 1 - b_allele_freq ;
		}
		
		result$a_allele_freq[i] = a_allele_freq
		result$b_allele_freq[i] = b_allele_freq

		result$missing_proportion[i] = ( nrow(G) - sum( G ) ) / nrow( G )
	}
	return( result ) ;
}



simulate.SNPs <- function( NSNPs, N_samples, case_control, sex, chromosome, model = "normal" ) {
	frequencies = runif( n = NSNPs )
	odds_ratios = exp( rnorm( n = NSNPs, mean = 0, sd = 0.5 ))
	null_SNPs = sample( x = c( TRUE, FALSE ), prob = c( 0.8, 0.2 ), replace = T, size = NSNPs )
	odds_ratios[ null_SNPs ] = 1

	model = rep( model, NSNPs )

	SNPs = data.frame(
		chromosome = chromosome,
		SNPID = sprintf( "SNP_%03d", 1:NSNPs ),
		rsid = sprintf( "rs%03d", 1:NSNPs ),
		position = 1:NSNPs,
		A_allele = 'A',
		B_allele = 'G',
		control_A_freq = 1 - frequencies,
		control_B_freq = frequencies,
		odds_ratio = odds_ratios,
		model = model
	)

	SNPs$case_B_freq = ( SNPs$control_B_freq / ( 1 - SNPs$control_B_freq ) ) * SNPs$odds_ratio
	SNPs$case_B_freq = SNPs$case_B_freq / ( 1 + SNPs$case_B_freq )
	SNPs$case_A_freq = 1 - SNPs$case_B_freq

	SNPs$control_AA_freq = (1 - SNPs$control_B_freq )^2
	SNPs$control_AB_freq = 2 * SNPs$control_B_freq * (1 - SNPs$control_B_freq )
	SNPs$control_BB_freq = SNPs$control_B_freq^2
	
	if( chromosome == "0X" && model != "single-active-copy" ) {
		diploid.OR = sqrt( SNPs$odds_ratio )
	} else {
		diploid.OR = SNPs$odds_ratio
	}
	SNPs$case_AA_freq = 1 / ( 1 + ( diploid.OR * SNPs$control_AB_freq / SNPs$control_AA_freq ) + ( diploid.OR^2 * SNPs$control_BB_freq / SNPs$control_AA_freq ) ) ;
	SNPs$case_AB_freq = SNPs$case_AA_freq * diploid.OR * SNPs$control_AB_freq / SNPs$control_AA_freq ;
	SNPs$case_BB_freq = SNPs$case_AA_freq * diploid.OR^2 * SNPs$control_BB_freq / SNPs$control_AA_freq ;

	genotypes = matrix( data = NA, nrow = nrow( SNPs ), ncol = 3 * N_samples )
	for( i in 1:nrow( SNPs )) {
		cat( i, "of", nrow(SNPs), "...\n" )
		genotypes[i,] = as.numeric(
			t(
				simulate.SNP.genotypes(
					case_control,
					sex,
					haploid_control_frequencies = c( SNPs$control_A_freq[i], SNPs$control_B_freq[i] ),
					haploid_case_frequencies = c( SNPs$case_A_freq[i], SNPs$case_B_freq[i] ),
					diploid_control_frequencies = SNPs[i, c( "control_AA_freq", "control_AB_freq", "control_BB_freq" ) ],
					diploid_case_frequencies = SNPs[i, c( "case_AA_freq", "case_AB_freq", "case_BB_freq" ) ],
					chromosome = SNPs$chromosome[i],
					model = model,
					debug = T
				)
			)
		) ;
	}
	return( list( id_data = SNPs, genotypes = genotypes ) ) ;
}

simulate.SNP.genotypes <- function(
	case_control, 					# a vector of zeroes and ones
	sex,          					# vector with "male"s and "female"s
	haploid_control_frequencies, 	# frequency of
	haploid_case_frequencies,
	diploid_control_frequencies,
	diploid_case_frequencies,
	chromosome,
	model = "normal", # or "single-active-copy" for single active copy model on X chromosome
	debug = FALSE
) {
	stopifnot( length( case_control ) == length( sex )) ;
	N_samples = length( case_control ) ;

	females = ( sex == "female" ) ;
	males = ( sex == "male" ) ;

	if( chromosome != "0X" ) {
		diploids = rep( TRUE, N_samples )
		haploids = rep( FALSE, N_samples )
	} else {
		diploids = females ;
		haploids = males ;
	}

	sample.haploid.control.genotypes <- function( i ) {
		if( is.na( case_control[i] ) ) {
			return( NA ) ;
		}
		return( which( rmultinom( n = 1, size = 1, prob = haploid_control_frequencies ) != 0 ) - 1 ) ;
	}

	sample.haploid.genotypes <- function( i ) {
		if( is.na( case_control[i] )) {
			return( NA )  ;
		} else if( case_control[i] == 1 ) {
			result = which( rmultinom( n = 1, size = 1, prob = haploid_case_frequencies ) != 0 ) ;
		} else {
			result = which( rmultinom( n = 1, size = 1, prob = haploid_control_frequencies ) != 0 ) ;
		}
		return( result - 1 ) ;
	}

	sample.diploid.genotypes <- function( i ) {
		if( is.na( case_control[i] )) {
			return( NA )  ;
		} else if( case_control[i] == 1 ) {
			result = which( rmultinom( n = 1, size = 1, prob = diploid_case_frequencies ) != 0 ) ;
		} else {
			result = which( rmultinom( n = 1, size = 1, prob = diploid_control_frequencies ) != 0 )
		}
		return( result - 1 ) ;
	}

	haploid.genotypes = sapply( 1:N_samples, sample.haploid.genotypes ) ;

	if( chromosome == "0X" && model == "single-active-copy" ) {
		diploid.genotypes = sapply( 1:N_samples, sample.haploid.genotypes ) ;
		diploid.genotypes = diploid.genotypes + sapply( 1:N_samples, sample.haploid.control.genotypes ) ;
	} else {
		diploid.genotypes = sapply( 1:N_samples, sample.diploid.genotypes ) ;
	}
	
	result = matrix( data = 0, ncol = 3, nrow = N_samples ) ;

	for( i in 1:N_samples ) {
		if( haploids[i] ) {
			result[ i, haploid.genotypes[i] + 1 ] = 1 ;
		}
		else if( diploids[i] ) {
			result[ i, diploid.genotypes[i] + 1 ] = 1 ;
		}
	}

	return( result )
}


# Now perform the simulation
save.sample.file( cohort[1:nrow(cohort1), ], filename = "cohort1.sample", column.types = c( "0", "0", "0", "D", "D", "C", "C", "P", "P", "B", "B", "D" ) )
save.sample.file( cohort[(nrow(cohort1)+1):nrow(cohort), ], filename = "cohort2.sample", column.types = c( "0", "0", "0", "D", "D", "C", "C", "P", "P", "B", "B", "D" ) )

case_control = cohort$bin2
sex = cohort$sex

N_samples = nrow( cohort )

snps = simulate.SNPs( 1000, nrow( cohort ), case_control, sex, chromosome = "01" )
write.table( cbind( snps$id_data[, 1:6], snps$genotypes[, 1:((N_samples/2)*3)] ), file = "cohort1_01.gen", col.names = F, row.names = F, quote = F )
write.table( cbind( snps$id_data[, 1:6], snps$genotypes[, (((N_samples/2)*3)+1):ncol( snps$genotypes )] ), file = "cohort2_01.gen", col.names = F, row.names = F, quote = F )
write.table( snps$id_data, file = "chromosome_01_snps.txt", col.names = T, row.names = F, quote = F )

snps = simulate.SNPs( 1000, nrow( cohort ), case_control, sex, "0X", model = "normal" )
write.table( cbind( snps$id_data[, 1:6], snps$genotypes[, 1:((N_samples/2)*3)] ), file = "cohort1_0X.gen", col.names = F, row.names = F, quote = F )
write.table( cbind( snps$id_data[, 1:6], snps$genotypes[, (((N_samples/2)*3)+1):ncol( snps$genotypes )] ), file = "cohort2_0X.gen", col.names = F, row.names = F, quote = F )
write.table( snps$id_data, file = "chromosome_0X_snps.txt", col.names = T, row.names = F, quote = F )

for( i in 1:N_samples ) {
	if( sex[i] == 'male' ) {
		snps$genotypes[, (3*i) ] = snps$genotypes[, (3*i)-1 ]
		snps$genotypes[, (3*i) - 1 ] = 0
	}
}
write.table( cbind( snps$id_data[, 1:6], snps$genotypes[, 1:((N_samples/2)*3)] ), file = "cohort1_0X_for_snptest.gen", col.names = F, row.names = F, quote = F )
write.table( cbind( snps$id_data[, 1:6], snps$genotypes[, (((N_samples/2)*3)+1):ncol( snps$genotypes )] ), file = "cohort2_0X_for_snptest.gen", col.names = F, row.names = F, quote = F )

snps = simulate.SNPs( 1000, nrow( cohort ), case_control, sex, "0X", model = "single-active-copy" )
write.table( cbind( snps$id_data[, 1:6], snps$genotypes[, 1:((N_samples/2)*3)] ), file = "cohort1_0X_single_active_copy.gen", col.names = F, row.names = F, quote = F )
write.table( cbind( snps$id_data[, 1:6], snps$genotypes[, (((N_samples/2)*3)+1):ncol( snps$genotypes )] ), file = "cohort2_0X_single-active-copy.gen", col.names = F, row.names = F, quote = F )
write.table( snps$id_data, file = "chromosome_0X_single-active-copy-snps.txt", col.names = T, row.names = F, quote = F )
