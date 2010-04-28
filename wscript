import os.path
import glob
import UnitTest
import Options

srcdir="."
APPNAME = "gen-tools"
VERSION = "1.2"

def set_options( opt ):
	opt.tool_options( 'compiler_cxx' )
	opt.tool_options( 'boost' )
	opt.add_option( "--staticbuild", action='store_true', default=False, help='Create statically-linked executables if possible.')

#-----------------------------------
# CONFIGURE
#-----------------------------------

def configure( conf ):
	conf.check_tool( 'compiler_cxx ')
	
	platform_specific_configure( conf )
	check_for_3rd_party_components( conf )
	misc_configure( conf )

	create_variant( conf, 'release' )
	configure_variant( conf, 'default', get_cxx_flags( 'default' ), get_ld_flags( 'default' ))
	configure_variant( conf, 'release', get_cxx_flags( 'release' ), get_ld_flags( 'release' ))

def create_variant( conf, variant_name ):
	variant = conf.env.copy()
	conf.set_env_name( variant_name, variant )
	variant.set_variant( variant_name )

def configure_variant( conf, variant_name, cxxflags, ldflags ):
	conf.setenv( variant_name )
	conf.env[ 'CXXFLAGS' ] = cxxflags
	conf.env[ 'LINKFLAGS' ] = ldflags
	conf.write_config_header( 'config.hpp' )
	conf.write_config_header( 'genfile/config.hpp' )

def check_for_3rd_party_components( conf ):
	check_for_boost_components( conf )
	check_for_zlib( conf )

def check_for_boost_components( conf ):
	conf.check_tool( 'boost' )
	if check_for_boost_headers( conf, '1.36.1' ):
		check_for_boost_lib( conf, 'iostreams', min_version='1.36', uselib="BOOST_IOSTREAMS" )
		check_for_boost_lib( conf, 'filesystem', min_version='1.36', uselib="BOOST_FILESYSTEM" )
		check_for_boost_lib( conf, 'system', min_version='1.36', uselib="BOOST_SYSTEM" )

def check_for_boost_headers( conf, min_version ):
	if conf.check_boost( min_version='1.36.1' ):
		conf.define( 'HAVE_BOOST_TIMER', 1 )
		conf.define( 'HAVE_BOOST_MATH', 1 )
		conf.define( 'HAVE_BOOST_FUNCTION', 1 )
		return True
	return False

def check_for_boost_lib( conf, lib, min_version, uselib ):
	if Options.options.staticbuild:
		static_selector = 'onlystatic'
	else:
		static_selector = 'nostatic'
	if conf.check_boost( lib = lib, min_version = min_version, static=static_selector, uselib = uselib ):
		conf.define( 'HAVE_' + uselib, 1 )

def check_for_zlib( conf ):
	if Options.options.staticbuild:
		if conf.check_cxx( staticlib='z', uselib_store='ZLIB' ):
			conf.define( 'HAVE_ZLIB', 1 )
	else:
		if conf.check_cxx( lib='z', uselib_store='ZLIB' ):
			conf.define( 'HAVE_ZLIB', 1 )

def platform_specific_configure( conf ):
	import platform
	if( platform.system() == 'Darwin' ):
		pass

def misc_configure( conf ) :
	conf.define ( 'GENFILE_USE_FAST_PARSE_METHODS', 1 )

def get_cxx_flags( variant_name ):
	cxxflags = ['-Wall', '-pedantic']
	if variant_name == 'default':
		cxxflags.extend( ['-g', '-p' ])
	elif variant_name == 'release':
		cxxflags.extend( [ '-O3' ])
	return cxxflags

def get_ld_flags( variant_name ):
	ldflags = []
	if Options.options.staticbuild:
		ldflags += [ '-static' ]
	return ldflags

#-----------------------------------
# BUILD
#-----------------------------------

def build( bld ):
	import Options
	
	bld.add_subdirs( 'genfile' )
	bld.add_subdirs( 'statfile' )
	
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
			'src/CmdLineOptionProcessor.cpp',
			'src/InformationStatistic.cpp',
			'src/SNPIDMatchesCondition.cpp'
		],
		includes='./include ./genfile/include',
		uselib = 'BOOST BOOST_IOSTREAMS ZLIB BOOST_MATH BOOST_FILESYSTEM BOOST_SYSTEM'
	)

	#---------------------
	# programs
	#---------------------
	create_app( bld, name='qctool', uselib_local = 'gen-tools-optionprocessor gen-tools-lib gen-tools-string gen-tools-exception genfile statfile' )
	create_app( bld, name='gen-convert', uselib_local = 'gen-tools-optionprocessor gen-tools-string gen-tools-exception gen-tools-lib genfile' )
	create_app( bld, name='gen-compare', uselib_local = 'gen-tools-optionprocessor gen-tools-string gen-tools-exception gen-tools-lib genfile' )
	create_app( bld, name='gen-grep', uselib_local = 'gen-tools-optionprocessor gen-tools-string gen-tools-exception gen-tools-lib genfile' )

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
	# misc tests...
	create_test( bld, 'test_log_of_gamma' )
	create_test( bld, 'test_log_of_factorial' )
	create_test( bld, 'test_genrow' )
	create_test( bld, 'test_genrow_io' )
	create_test( bld, 'test_row_conditions' )
	create_test( bld, 'test_wildcard' )
	# Statistic tests...
	create_test( bld, 'test_hwe' )
	create_test( bld, 'test_maximum_likelihood_statistics' )
	create_test( bld, 'test_hwe_against_SNPHWE' )
	create_test( bld, 'test_maf' )
	create_test( bld, 'test_missing' )
	create_test( bld, 'test_heterozygosity' )
	create_test( bld, 'test_alleles' )
	create_test( bld, 'test_statistic_arithmetic' )
	create_test( bld, 'test_information' )

def create_app( bld, name, uselib = '', uselib_local = '' ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'apps/' + name + '.cpp' ],
		includes='./include ./genfile/include ./statfile/include',
		uselib_local = uselib_local,
		uselib = uselib
	)

def create_test( bld, name ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'test/' + name + '.cpp' ],
		uselib_local = 'gen-tools-optionprocessor gen-tools-lib gen-tools-string gen-tools-exception genfile',
		includes='./include ./genfile/include',
		unit_test=1
	)

def create_benchmark( bld, name ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'benchmarks/' + name + '.cpp' ],
		uselib_local = 'gen-tools-optionprocessor gen-tools-lib gen-tools-string gen-tools-exception genfile',
		includes='./include ./genfile/include'
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
