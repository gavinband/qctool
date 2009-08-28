import os.path
import glob
import UnitTest

srcdir="."

def get_version():
	return "1.0_beta4"

def set_options( opt ):
	opt.tool_options( 'compiler_cxx' )

#-----------------------------------
# CONFIGURE
#-----------------------------------

def configure( conf ):
	conf.check_tool( 'compiler_cxx ')
	
	platform_specific_configure( conf )
	check_for_3rd_party_components( conf )
	misc_configure( conf )

	create_variant( conf, 'release' )
	configure_variant( conf, 'default', ['-g', '-p', '-Wall'] )
	configure_variant( conf, 'release', ['-Wall', '-O3'] )

def create_variant( conf, variant_name ):
	variant = conf.env.copy()
	conf.set_env_name( variant_name, variant )
	variant.set_variant( variant_name )

def configure_variant( conf, variant_name, cxxflags ):
	conf.setenv( variant_name )
	conf.env[ 'CXXFLAGS' ] = cxxflags
	conf.write_config_header( 'config.hpp' )
	conf.write_config_header( 'genfile/config.hpp' )

def check_for_3rd_party_components( conf ):
	check_for_boost_components( conf )
	check_for_zlib( conf )

def check_for_boost_components( conf ):
	conf.check_tool( 'boost' )

	if conf.check_boost( min_version='1.36.1' ):
		conf.define( 'HAVE_BOOST_TIMER', 1 )
		conf.define( 'HAVE_BOOST_MATH', 1 )
		conf.define( 'HAVE_BOOST_FUNCTION', 1 )
	if conf.check_boost( lib='iostreams', min_version='1.36', uselib="BOOST_IOSTREAMS" ):
		conf.define( 'HAVE_BOOST_IOSTREAMS', 1 )
	if conf.check_boost( lib='filesystem', min_version='1.36', uselib="BOOST_FILESYSTEM" ):
		conf.define( 'HAVE_BOOST_FILESYSTEM', 1 )
	if conf.check_boost( lib='system', min_version='1.36', uselib="BOOST_SYSTEM" ):
		conf.define( 'HAVE_BOOST_SYSTEM', 1 )
	if conf.check_boost( lib='random', min_version='1.36', uselib='BOOST_RANDOM' ):
		conf.define( 'HAVE_BOOST_RANDOM', 1 )

def check_for_zlib( conf ):
	if conf.check_cxx( lib='z', uselib_store='ZLIB' ):
		conf.define( 'HAVE_ZLIB', 1 )

def platform_specific_configure( conf ):
	import platform
	if( platform.system() == 'Darwin' ):
		pass

def misc_configure( conf ) :
	conf.define ( 'GENFILE_USE_FAST_PARSE_METHODS', 1 )


#-----------------------------------
# BUILD
#-----------------------------------

def build( bld ):
	import Options
	
	bld.add_subdirs( 'genfile' )
	
	#---------------------
	# libs
	#---------------------
	bld.new_task_gen(
		features = 'cxx cstaticlib',
		target = 'gen-tools-exception',
		source = [  'src/GToolException.cpp' ],
		includes='./include'
	)

	bld.new_task_gen(
		features = 'cxx cstaticlib',
		target = 'gen-tools-string',
		source = [  'src/string_utils.cpp' ],
		includes='./include'
	)

	bld.new_task_gen(
		features = 'cxx cstaticlib',
		target = 'gen-tools-optionprocessor',
		source = [  'src/OptionProcessor.cpp',
		 			'src/OptionDefinition.cpp'
		],
		includes='./include',
		uselib = 'gen-tools-exception gen-tools-string'
	)

	bld.new_task_gen(
		features = 'cxx cstaticlib',
		target = 'gen-tools-lib',
		source = [  
			'src/AlleleProportions.cpp',
			'src/Condition.cpp',
			'src/FileUtil.cpp',
			'src/wildcard.cpp',
			'src/GToolException.cpp',
			'src/GenRow.cpp',
			'src/GenRowIO.cpp',
			'src/ExternalStorageGenRow.cpp',
			'src/GenRowStatistics.cpp',
			'src/GenotypeAssayBasicStatistics.cpp',
			'src/GenotypeAssayStatisticArithmetic.cpp',
			'src/GenotypeAssayStatisticFactory.cpp',
			'src/GenotypeAssayStatistics.cpp',
			'src/GenotypeProportions.cpp',
			'src/HardyWeinbergExactTestStatistic.cpp',
			'src/LikelihoodRatioTestStatistic.cpp',
			'src/RowCondition.cpp',
			'src/SNPHWE.cpp',
			'src/SNPInListCondition.cpp',
			'src/SampleRow.cpp',
			'src/SampleRowStatistics.cpp',
			'src/SimpleGenotypeAssayStatistics.cpp',
			'src/SampleInListCondition.cpp',
			'src/Whitespace.cpp',
			'src/distributions.cpp',
			'src/gamma.cpp',
			'src/parse_utils.cpp',
			'src/string_to_value_map.cpp',
			'src/FileBackupCreator.cpp',
			'src/InputToOutputFilenameMapper.cpp',
			'src/OstreamTee.cpp',
			'src/CmdLineOptionProcessor.cpp'
		],
		includes='./include ./genfile/include',
		uselib = 'BOOST BOOST_IOSTREAMS BOOST_MATH BOOST_FILESYSTEM BOOST_SYSTEM'
	)

	#---------------------
	# programs
	#---------------------
	create_app( bld, name='qc-tool', uselib_local = 'gen-tools-lib gen-tools-string gen-tools-exception gen-tools-optionprocessor genfile', uselib = 'RLIB' )
	create_app( bld, name='gen-convert', uselib_local = 'gen-tools-lib gen-tools-string gen-tools-exception gen-tools-optionprocessor genfile' )
	create_app( bld, name='gen-compare', uselib_local = 'gen-tools-lib gen-tools-string gen-tools-exception gen-tools-optionprocessor genfile' )
	create_app( bld, name='gen-case-control-test', uselib_local = 'gen-tools-lib gen-tools-string gen-tools-exception gen-tools-optionprocessor genfile' )
	create_app( bld, name='generate-random-permutations-of-0-1-vector', uselib_local = 'gen-tools-string gen-tools-exception gen-tools-optionprocessor', uselib = 'BOOST BOOST_RANDOM' )

	#---------------------
	# benchmarks
	#---------------------
	create_benchmark( bld, 'benchmark-statistics' )

	# Build release variants as well.
	for obj in [] + bld.all_task_gen:
	    new_obj = obj.clone('release')

	#---------------------
	# tests
	#---------------------
	create_test( bld, 'test_log_of_gamma' )
	create_test( bld, 'test_log_of_factorial' )
	create_test( bld, 'test_genrow' )
	create_test( bld, 'test_genrow_io' )
	create_test( bld, 'test_genbin_snp_format' )
	create_test( bld, 'test_hardy_weinberg_exact_test_statistic' )
	create_test( bld, 'test_maximum_likelihood_statistics' )
	create_test( bld, 'test_hardy_weinberg_exact_test_statistic_against_SNPHWE' )
	create_test( bld, 'test_maf' )
	create_test( bld, 'test_missing' )
	create_test( bld, 'test_heterozygosity' )
	create_test( bld, 'test_alleles' )
	create_test( bld, 'test_statistic_arithmetic' )
	create_test( bld, 'test_row_conditions' )
	create_test( bld, 'test_wildcard' )


def create_app( bld, name, uselib = '', uselib_local = '' ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'apps/' + name + '.cpp' ],
		includes='./include ./genfile/include',
		uselib_local = uselib_local,
		uselib = uselib
	)

def create_test( bld, name ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'test/' + name + '.cpp' ],
		uselib_local = 'gen-tools-lib gen-tools-string gen-tools-exception gen-tools-optionprocessor genfile',
		includes='./include ./genfile/include',
		unit_test=1
	)

def create_benchmark( bld, name ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'benchmarks/' + name + '.cpp' ],
		uselib_local = 'gen-tools-lib gen-tools-string gen-tools-exception gen-tools-optionprocessor genfile',
		includes='./include ./genfile/include',
		unit_test=1
	)

#-----------------------------------
# CHECK
#-----------------------------------

def check(context):
	# Unit tests are run when "check" target is used
	ut = UnitTest.unit_test()
	ut.change_to_testfile_dir = True
	ut.run()
	ut.print_results()
