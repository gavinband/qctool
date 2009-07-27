import os.path
import glob
import UnitTest

srcdir="."


def set_options( opt ):
	opt.tool_options( 'compiler_cxx' )

#-----------------------------------
# CONFIGURE
#-----------------------------------

def configure( conf ):
	conf.check_tool( 'compiler_cxx ')

	platform_specific_configure( conf )
	check_for_boost_components( conf )
	check_for_zlib( conf )
	misc_configure( conf )

	create_variant( conf, 'release' )
	create_default_variant( conf )
	create_release_variant( conf )

def create_variant( conf, variant_name ):
	variant = conf.env.copy()
	conf.set_env_name( variant_name, variant )
	variant.set_variant( variant_name )

def create_default_variant( conf ):
	conf.setenv( 'default')
	conf.env.append_value( 'CXXFLAGS', [ '-g', '-Wall', '-p' ] )
	conf.write_config_header( 'config.hpp' )

def create_release_variant( conf ):
	conf.setenv( 'release' )
	conf.env.append_value( 'CXXFLAGS', '-Wall' ) 
	conf.write_config_header( 'config.hpp' )

def check_for_3rd_part_components( conf ):
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
	if conf.check_cxx( lib = 'z', uselib_store='ZLIB' ):
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
		target = 'gtool-exception',
		source = [  'src/GToolException.cpp' ],
		includes='./include'
	)

	bld.new_task_gen(
		features = 'cxx cstaticlib',
		target = 'gtool-optionprocessor',
		source = [  'src/OptionProcessor.cpp',
		 			'src/OptionDefinition.cpp'
		],
		includes='./include',
		uselib = 'gtool-exception'
	)

	bld.new_task_gen(
		features = 'cxx cstaticlib',
		target = 'gtool-lib',
		source = [  
			'src/AlleleProportions.cpp',
			'src/Condition.cpp',
			'src/FileUtil.cpp',
			'src/wildcard.cpp',
			'src/GToolException.cpp',
			'src/GenRow.cpp',
			'src/GenRowIO.cpp',
			'src/ExternalStorageGenRow.cpp',
			'src/GenRowSource.cpp',
			'src/GenRowSink.cpp',
			'src/GenRowStatistics.cpp',
			'src/GenotypeAssayBasicStatistics.cpp',
			'src/GenotypeAssayStatisticArithmetic.cpp',
			'src/GenotypeAssayStatisticFactory.cpp',
			'src/GenotypeAssayStatisticInRange.cpp',
			'src/GenotypeAssayStatistics.cpp',
			'src/GenotypeProportions.cpp',
			'src/HardyWeinbergExactTestStatistic.cpp',
			'src/LikelihoodRatioTestStatistic.cpp',
			'src/OptionProcessor.cpp',
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
			'src/string_utils.cpp',
			'src/string_to_value_map.cpp'
		],
		includes='./include ./genfile/include',
		uselib = 'BOOST BOOST_IOSTREAMS BOOST_MATH BOOST_FILESYSTEM BOOST_SYSTEM'
	)

	#---------------------
	# programs
	#---------------------
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = 'gen-select',
		source = [  'apps/gen-select.cpp' ],
		includes='./include ./genfile/include',
		uselib_local = 'gtool-lib gtool-exception gtool-optionprocessor genfile'
	)

	bld.new_task_gen(
		features = 'cxx cprogram',
		target = 'gen-convert',
		source = [  'apps/gen-convert.cpp' ],
		includes='./include ./genfile/include',
		uselib_local = 'gtool-lib gtool-exception gtool-optionprocessor genfile'
	)

	bld.new_task_gen(
		features = 'cxx cprogram',
		target = 'gen-compare',
		source = [  'apps/gen-compare.cpp' ],
		includes='./include ./genfile/include',
		uselib_local = 'gtool-lib gtool-exception gtool-optionprocessor genfile'
	)

	bld.new_task_gen(
		features = 'cxx cprogram',
		target = 'gen-case-control-test',
		source = [  'apps/gen-case-control-test.cpp' ],
		includes='./include ./genfile/include',
		uselib_local = 'gtool-lib gtool-exception gtool-optionprocessor genfile'
	)

	bld.new_task_gen(
		features = 'cxx cprogram',
		target = 'generate-random-permutations-of-0-1-vector',
		source = [  'apps/generate-random-permutations-of-0-1-vector.cpp' ],
		includes='./include ./genfile/include',
		uselib_local = 'gtool-exception gtool-optionprocessor',
		uselib = 'BOOST BOOST_RANDOM'
	)

	#---------------------
	# benchmarks
	#---------------------
	create_benchmark( bld, 'benchmark-statistics' )
	create_benchmark( bld, 'benchmark-io' )

	# Build release variants as well.
	for obj in [] + bld.all_task_gen:
	    new_obj = obj.clone('release')

	#---------------------
	# tests
	#---------------------
	create_test( bld, 'test_log_of_gamma' )
	create_test( bld, 'test_log_of_factorial' )
	create_test( bld, 'test_genrow_io' )
	create_test( bld, 'test_genbin_snp_format' )
	create_test( bld, 'test_hardy_weinberg_exact_test_statistic' )
	create_test( bld, 'test_maximum_likelihood_statistics' )
	create_test( bld, 'test_hardy_weinberg_exact_test_statistic_against_SNPHWE' )
	create_test( bld, 'test_simple_statistics' )
	create_test( bld, 'test_statistic_arithmetic' )
	create_test( bld, 'test_fileutil' )


def create_test( bld, name ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'test/' + name + '.cpp' ],
		uselib_local = 'gtool-lib gtool-exception gtool-optionprocessor genfile',
		includes='./include ./genfile/include',
		unit_test=1
	)

def create_benchmark( bld, name ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'benchmarks/' + name + '.cpp' ],
		uselib_local = 'gtool-lib gtool-exception gtool-optionprocessor genfile',
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
