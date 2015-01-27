library( hapdb )
db = open.hapdb( "~/Projects/Software/snptest/example/cohort1.sqlite" )
system.time( { D2 = load.genotypes( db, chromosome = "01", range = c( 0, 1000), method = "cpp", compute.dosage = TRUE ) })
system.time( { D1 = load.genotypes( db, chromosome = "01", range = c( 0, 1000), method = "R", compute.dosage = TRUE ) })

options(width=200)
D1$data[1:10,1:15]
D2$data[1:10,1:15]

D1$dosage[1:10,1:15]
D2$dosage[1:10,1:15]

stopifnot( length( which( D1$data != D2$data )) == 0 )
stopifnot( length( which( D1$dosage != D2$dosage )) == 0 )
