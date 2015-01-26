#!/bin/bash
ORIGINAL=`pwd`

OUTCOME=test_outcomes/7d74cae1bc07736decac3c4baedee607
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example_#.gen -og ./example.bgen
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME};

OUTCOME=test_outcomes/e2895a2ce9da2521a90a0c2895a07913
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example_#.gen -og ./example.vcf 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME};

OUTCOME=test_outcomes/b7145d3d59dec793f4c736a1326d645b
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example_#.gen -og ./example_#.bgen 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME};

OUTCOME=test_outcomes/4963d28ec976cd2defc284bc67230c40
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example_01.gen -og ./example_01.bgen -assume-chromosome 01 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/5d18c029e59537ccf913ad61462355f3
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example_01.gen.gz -og example_01.bgen -assume-chromosome 01
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/10717346bcd1300044d6f1e059b04a48
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.vcf -og ./converted.bgen 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/c9303cf8c3848570ae21979089abf2a4
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.vcf.gz -og converted.bgen 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/ee429ab306285b26e984149c559abc74
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.4 -g ${ORIGINAL}/setup/example_bgzipped.vcf.gz -og converted.bgen
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/f1e75291f44c79b76f7a225a98c3971f
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.vcf -vcf-genotype-field GP -og ./converted.bgen 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/4e58b1968cd59121e5e188f5ce3a1122
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -sort -og ./sorted.bgen 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/beb81eda3016d5a7e295fc48a0547cd4
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -og ./subsetted.gen -incl-rsids ${ORIGINAL}/setup/rsids_to_include.txt
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/2ec7633648c7d0cfba37c3022fe7944c
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -og ./subsetted.gen -excl-range 06:15000-18000
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/016b47e1e61b919262b78a19077d491e
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -og ./subsetted.gen -snp-missing-rate 0.05 -maf 0.01 1 -info 0.9 1 -hwe 20 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/7afeb9f641787813e71b33c68a1a43c7
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -write-snp-excl-list ./snp_exclusions.txt -snp-missing-rate 0.05 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/dfb2160f47f237e68d7a419f92375db4
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -s ${ORIGINAL}/setup/example.sample -sample-stats ./example.sample-stats -os ./example.sample 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/4ad6e5d7a116f0020f04557227391049
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -snp-stats ./example.snp-stats 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/55acb5cfea24aa16b9339eb91e88ba2b
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -s ${ORIGINAL}/setup/example.sample -snp-stats example.snp-stats -excl-samples ${ORIGINAL}/setup/samples_to_exclude.txt 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/89cd8cb1946efb9ec1aff5ca3ab29c2b
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -s ${ORIGINAL}/setup/example.sample -g ${ORIGINAL}/setup/second_cohort.bgen -s ${ORIGINAL}/setup/second_cohort.sample -og ./joined.gen -os ./joined.sample 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/8ef03f6764face5b96053304632369f1
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -s ${ORIGINAL}/setup/example.sample -g ${ORIGINAL}/setup/second_cohort.bgen -s ${ORIGINAL}/setup/second_cohort.sample -og ./joined.gen -os ./joined.sample -snp-match-fields position,alleles 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/a24a80ab653d479d7320588b713572ee
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -s ${ORIGINAL}/setup/example.sample -g ${ORIGINAL}/setup/second_cohort.bgen -s ${ORIGINAL}/setup/second_cohort.sample -og ./joined.gen -os ./joined.sample -snp-match-fields position,alleles -match-alleles-to-cohort1 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/d64e4342281c2124590073cd88610eb4
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example_#.gen -s ${ORIGINAL}/setup/example.sample -os ./example.sample -condition-on rs1234 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/f6315fbf2c41fb862202c3ca6bdcd33d
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example_#.gen -s ${ORIGINAL}/setup/example.sample -os ./example.sample -condition-on pos~03:10001 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/0fc08af9c3233157dbb613a96ee495ba
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -s ${ORIGINAL}/setup/example.sample -os ./example.sample -quantile-normalise pheno1 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}

OUTCOME=test_outcomes/08779fe4397fd269e5e674acd0429928
mkdir -p ${OUTCOME}; cd ${OUTCOME}
#tar ${ORIGINAL}/setup.tgz
qctool-1.3 -g ${ORIGINAL}/setup/example.bgen -s ${ORIGINAL}/setup/example.sample -os ./example.sample -quantile-normalise pheno1,pheno2 
cd ${ORIGINAL}; chmod -R ugo-w ${OUTCOME}