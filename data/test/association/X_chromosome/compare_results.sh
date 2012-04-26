#!/bin/bash

rm test.qcdb
mkdir -p images

qctool-dev \
   -g cohort1_0X.gen \
   -s cohort1.sample \
   -g cohort2_0X.gen \
   -s cohort2.sample \
   -snp-stats test.qcdb \
   -test bin2 \
   -analysis-name "full"

qctool-dev \
  -g cohort1_0X_autosomalised.gen \
  -s cohort1.sample \
  -g cohort2_0X_autosomalised.gen \
  -s cohort2.sample \
  -snp-stats test.qcdb \
  -test bin2 \
  -analysis-name "autosomalised"

R --vanilla << HERE_DOC

library( GWAS )
library( RSQLite )

expected = read.table( "X_chromosome_SNP_information.txt", hea=T, as.is=T )
db = dbConnect( dbDriver( "SQLite" ), dbname = "test.qcdb" )

dbGetQuery( db, "SELECT * FROM Entity" )

results = dbGetQuery( db, "SELECT * FROM SummaryDataView WHERE name == 'full'" )

summary( expected )
summary( results )

beta_1 = which( results[, "variable" ] == "beta_1" )
p_value = which( results[, "variable" ] == "p_value" )
alleleB_frequency = which( results[, "variable" ] == "alleleB_frequency" )

png( file = "images/results_versus_expected_with_X_inactivation.png", width = 1200, height = 1200 )
par( mfrow = c(2,2) )
col = rep( "black", nrow( expected ))
col[ which( expected[, "odds_ratio" ] != 1 )] = "red"

betas = results[ beta_1, c( "rsid", "value" ) ]
p_values =  results[ p_value, c( "rsid", "value" ) ]
freqs = results[ alleleB_frequency, c( "rsid", "value" ) ]

M = match( expected[, "rsid" ], betas[, "rsid" ] )
plot( log( expected[, "odds_ratio" ] ), betas[ M, "value" ], xlab = "Expected log-odds ratio", ylab = "Fitted log-odds ratio", col = col )

M = match( expected[, "rsid" ], p_values[, "rsid" ] )
plot( 1:nrow( expected ), -log10( p_values[ M, "value" ] ), xlab = "SNP index", ylab = "-log10 P-value", col = col )

M = match( expected[, "rsid" ], freqs[, "rsid" ] )
plot( expected[, "control_frequency" ], freqs[ M, "value" ], xlab = "Expected frequency in controls", ylab = "Fitted frequency in all samples", col = col )
dev.off()

HERE_DOC