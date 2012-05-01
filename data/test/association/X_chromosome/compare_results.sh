#!/bin/bash

rm test.qcdb
mkdir -p images

~/Projects/Software/qctool/build/release/qctool-dev \
   -g cohort1_0X.gen \
   -s cohort1.sample \
   -g cohort2_0X.gen \
   -s cohort2.sample \
   -snp-stats test.qcdb \
   -test bin2 \
   -analysis-name "full"


snptest_v2.3.0 \
    -data cohort1_0X_for_snptest.gen \
    cohort1.sample \
    cohort2_0X_for_snptest.gen \
    cohort2.sample \
    -frequentist 1 \
    -method ml \
    -o test.snptest \
    -pheno bin2

R --vanilla << HERE_DOC

library( GWAS )
library( RSQLite )

expected = read.table( "chromosome_0X_snps.txt", hea=T, as.is=T )
db = dbConnect( dbDriver( "SQLite" ), dbname = "test.qcdb" )

dbGetQuery( db, "SELECT * FROM Entity" )

results = dbGetQuery( db, "SELECT * FROM SummaryDataView WHERE name == 'full'" )

snptest = read.snptest( "test.snptest" )

summary( expected )
summary( results )

beta_1 = which( results[, "variable" ] == "beta_1" )
p_value = which( results[, "variable" ] == "p_value" )
alleleB_frequency = which( results[, "variable" ] == "alleleB_frequency" )

png( file = "images/results_versus_expected_with_X_inactivation.png", width = 1200, height = 900 )
par( mfrow = c(2,3) )

pch = rep( "o", nrow( expected ))
pch[ which( expected[, "control_B_freq" ] < 0.01 | expected[, "control_B_freq" ] > 0.99 ) ] = "+"

col = rep( "black", nrow( expected ))
col[ which( expected[, "odds_ratio" ] != 1 )] = "red"

betas = results[ beta_1, c( "rsid", "value" ) ]
p_values =  results[ p_value, c( "rsid", "value" ) ]
freqs = results[ alleleB_frequency, c( "rsid", "value" ) ]

M = match( expected[, "rsid" ], betas[, "rsid" ] )
plot( log( expected[, "odds_ratio" ] ), betas[ M, "value" ], xlab = "Expected log-odds ratio", ylab = "Fitted log-odds ratio", col = col, pch = pch )

M = match( expected[, "rsid" ], p_values[, "rsid" ] )
plot( 1:nrow( expected ), -log10( p_values[ M, "value" ] ), xlab = "SNP index", ylab = "-log10 P-value", col = col, pch = pch )

M = match( expected[, "rsid" ], freqs[, "rsid" ] )
plot( expected[, "control_B_freq" ], freqs[ M, "value" ], xlab = "Expected frequency in controls", ylab = "Fitted frequency in all samples", col = col, pch = pch )

M = match( snptest[, "rsid" ], betas[, "rsid" ] )
plot( snptest[, "bin2_frequentist_add_ml_beta_1" ], betas[ M, "value" ], xlab = "log-odds ratio fitted by SNPTEST", ylab = "log-odds ratio fitted by qctool", col = col, pch = pch )

M = match( snptest[, "rsid" ], p_values[, "rsid" ] )
plot( -log10( snptest[, "bin2_frequentist_add_ml_pvalue" ] ), -log10( p_values[ M, "value" ] ), xlab = "-log10 P-valuefitted by SNPTEST", ylab = "-log10 P-value fitted by qctool", col = col, pch = pch )
dev.off()

HERE_DOC
