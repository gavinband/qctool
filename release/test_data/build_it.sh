#!/bin/bash
cd /Users/gav/Projects/Software/qctool/release/test_data/test

rm -f *
qctool-1.1 -g ../setup/example_#.gen -og ./example.bgen 
cd test
tar -czf ../e7e9be61e6d008b0076987496b80c45f.tgz ./*

rm -f *
qctool-1.1 -g ../setup/example_#.gen -og ./example.vcf 
tar -czf ../7d74cae1bc07736decac3c4baedee607.tgz ./*

rm -f *
qctool-1.1 -g ../setup/example_#.gen -og ./example_#.bgen 
tar -czf ../b7145d3d59dec793f4c736a1326d645b.tgz ./*

rm -f *
qctool-1.1 -g ../setup/example_01.gen -og ./example_01.bgen -assume-chromosome 01 
tar -czf ../4963d28ec976cd2defc284bc67230c40.tgz ./*

rm -f *
qctool-1.1 -g %INPUT_DIR%/example_01.gen.gz -og %OUTPUT_DIR%/example_01.bgen -assume-chromosome 01
tar -czf ../5d18c029e59537ccf913ad61462355f3.tgz ./*

rm -f *
qctool-1.1 -g ../setup/example.vcf -og ./converted.bgen 
tar -czf ../10717346bcd1300044d6f1e059b04a48.tgz ./*

rm -f *
-g %INPUT_DIR%/example.vcf.gz -og %OUTPUT_DIR%/converted.bgen 
tar -czf ../c9303cf8c3848570ae21979089abf2a4.tgz ./*

rm -f *
qctool-1.1 -g %INPUT_DIR%/example.vcf.bgz -og %OUTPUT_DIR%/converted.bgen
tar -czf ../49bde70f20e3f7379868803d7e7a5410.tgz ./*

rm -f *
qctool-1.1 -g ../setup/example.vcf -vcf-genotype-field GP -og ./converted.bgen 
tar -czf ../f1e75291f44c79b76f7a225a98c3971f.tgz ./*

rm -f *
qctool-1.1 -g ../setup/example.bgen -sort -og ./sorted.bgen 
tar -czf ../fdae6ae2cca6af6fe9055e87b7b617bb.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -og ./subsetted.gen -incl-rsids ../setup/rsids_to_include.txt
tar -czf ../02b4e7b9d19134bff475e4ab47774f02.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -og ./subsetted.gen -excl-range 06:15000-18000
tar -czf ../2ec7633648c7d0cfba37c3022fe7944c.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -og ./subsetted.gen -snp-missing-rate 0.05 -maf 0.01 1 -info 0.9 1 -hwe 20 
tar -czf ../016b47e1e61b919262b78a19077d491e.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -write-snp-excl-list ./snp_exclusions.txt -snp-missing-rate 0.05 
tar -czf ../7afeb9f641787813e71b33c68a1a43c7.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -s ../setup/example.sample -sample-stats ./example.sample-stats -os ./example.sample 
tar -czf ../dfb2160f47f237e68d7a419f92375db4.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -snp-stats ./example.snp-stats 
tar -czf ../4ad6e5d7a116f0020f04557227391049.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -s ../setup/example.sample -snp-stats example.snp-stats -excl-samples ../setup/samples_to_exclude.txt 
tar -czf ../55acb5cfea24aa16b9339eb91e88ba2b.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -s ../setup/example.sample -g ../setup/second_cohort.bgen -s ../setup/second_cohort.sample -og ./joined.gen -os ./joined.sample 
tar -czf ../89cd8cb1946efb9ec1aff5ca3ab29c2b.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -s ../setup/example.sample -g ../setup/second_cohort.bgen -s ../setup/second_cohort.sample -og ./joined.gen -os ./joined.sample -snp-match-fields position,alleles 
tar -czf ../8ef03f6764face5b96053304632369f1.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -s ../setup/example.sample -g ../setup/second_cohort.bgen -s ../setup/second_cohort.sample -og ./joined.gen -os ./joined.sample -snp-match-fields position,alleles -match-alleles-to-cohort1 
tar -czf ../a24a80ab653d479d7320588b713572ee.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example_#.gen -s ../setup/example.sample -os ./example.sample -condition-on rs1234 
tar -czf ../d64e4342281c2124590073cd88610eb4.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example_#.gen -s ../setup/example.sample -os ./example.sample -condition-on pos~03:10001 
tar -czf ../f6315fbf2c41fb862202c3ca6bdcd33d.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -s ../setup/example.sample -os ./example.sample -quantile-normalise pheno1 
tar -czf ../0fc08af9c3233157dbb613a96ee495ba.tgz ./*


rm -f *
qctool-1.1 -g ../setup/example.bgen -s ../setup/example.sample -os ./example.sample -quantile-normalise pheno1,pheno2 
tar -czf ../08779fe4397fd269e5e674acd0429928.tgz ./*
