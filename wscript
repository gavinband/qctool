import os.path
import glob
#import UnitTest
import Options

srcdir="."
APPNAME = "qctool"
VERSION = "2.1-dev"

subdirs = [
	'genfile', 'statfile', 'string_utils', 'appcontext',
	'fputils', 'worker', 'integration',
	'3rd_party', 'components', 'qcdb', 'metro'
]

def set_options( opt ):
	opt.tool_options( 'compiler_cxx' )
	opt.tool_options( 'compiler_cc' )
	#opt.tool_options( 'boost' )
	opt.add_option( "--static", action='store_true', default=False, help='Create statically-linked executables if possible.')
	opt.add_option( "--all_targets", action='store_true', default=False, help='Create all targets, not just qctool.')

#-----------------------------------
# CONFIGURE
#-----------------------------------

def configure( conf ):
	print( "Using prefix\t\t\t\t :", conf.env[ 'PREFIX' ] )

	conf.check_tool( 'compiler_cxx')
	conf.check_tool( 'compiler_cc')

	cxxflags = conf.env[ 'CXXFLAGS' ] 
	linkflags = conf.env[ 'LINKFLAGS' ]
	if conf.check( cxxflags = '-std=c++11' ):
		cxxflags.append( '-std=c++11' )
		cxxflags.append( '-Wno-deprecated-declarations' )
	else:
		cxxflags.append( '-std=c++98' )

	# Disabling vectorisation checks as not always safe.
	# Use CXXFLAGS environment variable to turn these on
	# for flag in [ '-msse2', '-mavx', '-mssse3']:
	#	if conf.check( cxxflags = flag ):
	#		cxxflags.append( flag )
	
	platform_specific_configure( conf )
	check_for_3rd_party_components( conf )
	misc_configure( conf )

	import platform
	if Options.options.static and platform.system() != "Darwin":
		conf.env.SHLIB_MARKER='-Wl,-Bstatic'
	create_variant( conf, 'release' )
	configure_variant( conf, 'default', cxxflags, linkflags )
	configure_variant( conf, 'release', cxxflags, linkflags )

def check_for_3rd_party_components( conf ):
	check_for_zlib( conf )
	conf.define( 'HAVE_SQLITE3', 1 )
	conf.define( 'HAVE_EIGEN', 1 )
	conf.define( 'HAVE_ZSTD', 1 )
	if conf.check_cxx( lib = 'dl', uselib_store = 'DL' ):
		conf.define( 'HAVE_DL', 1 )
	if conf.check_cxx( lib = 'rt', uselib_store = 'RT' ):
		conf.define( 'HAVE_RT', 1 )
	if conf.check_cxx( lib = 'm', uselib_store = 'M' ):
		conf.define( 'HAVE_M', 1 )

	# Remove support for some things
	conf.define( 'HAVE_BZIP2', 0 )
	conf.define( 'HAVE_MGL', 0 )
	conf.define( 'HAVE_CAIRO', 0 )

	if Options.options.static and conf.check_cxx( staticlib = 'pthread', uselib_store = "PTHREAD" ):
		conf.define( 'HAVE_PTHREAD', 1 )
	elif conf.check_cxx( lib = 'pthread', uselib_store = "PTHREAD" ):
		conf.define( 'HAVE_PTHREAD', 1 )
	if conf.check_cc( lib = 'readline', uselib_store = 'READLINE' ):
		conf.define( "HAVE_READLINE", 1 )
	
	# Boost libs are now contained in the repo
	conf.define( "HAVE_BOOST_IOSTREAMS", 1 )
	conf.define( "HAVE_BOOST_FILESYSTEM", 1 )
	conf.define( "HAVE_BOOST_SYSTEM", 1 )
	conf.define( "HAVE_BOOST_THREAD", 1 )
	conf.define( "HAVE_BOOST_DATE_TIME", 1 )
	conf.define( "HAVE_BOOST_UNIT_TEST_FRAMEWORK", 1 )
	conf.define( "HAVE_BOOST_TIMER", 1 )
	conf.define( "HAVE_BOOST_REGEX", 1 )
	conf.define( "HAVE_BOOST_MATH", 1 )
	conf.define( "HAVE_BOOST_FUNCTION", 1 )
	conf.define( "HAVE_BOOST_SPIRIT", 1 )
	conf.define( "BOOST_FILESYSTEM_VERSION", 3 )

def check_for_zlib( conf ):
	if conf.check_cxx( staticlib='z', uselib_store='ZLIB' ):
		conf.define( 'HAVE_ZLIB', 1 )
	elif conf.check_cxx( lib='z', uselib_store='ZLIB' ):
		conf.define( 'HAVE_ZLIB', 1 )

def platform_specific_configure( conf ):
	import platform
	if platform.system() == 'Darwin':
		if conf.check_cxx( header_name='mach/mach_time.h', uselib_store = 'MACH_TIME' ):
			conf.define( 'HAVE_MACH_TIME', 1 )
		if conf.check_cxx(
			lib = 'cblas',
			fragment = '#include "cblas.h"\nint main() {}',
			cxxflags = '-I/System/Library/Frameworks/vecLib.framework/Headers',
			uselib_store = 'CBLAS'
		):
			conf.define( 'HAVE_CBLAS', 1 )
		if conf.check_cxx(
			header_name = 'clapack.h',
			lib = 'clapack',
			cxxflags = '-I/System/Library/Frameworks/vecLib.framework/Headers',
			uselib_store = 'CLAPACK'
		):
			conf.define( 'HAVE_CLAPACK', 1 )
			conf.define( 'HAVE_LAPACK', 1 )
	else:
		blasCode = """#include "cblas.h"
int main() {
	float v[10] ;
	float result = cblas_snrm2( 1, v, 1 ) ;
	return int(result) ;
}
"""
		if conf.check_cxx( lib = 'cblas', fragment = blasCode, uselib_store = 'CBLAS' ):
			conf.define( 'HAVE_CBLAS', 1 ) ;
		elif conf.check_cxx( lib = 'blas', fragment = blasCode, uselib_store = 'CBLAS' ):
			conf.define( 'HAVE_CBLAS', 1 ) ;
		elif conf.check_cxx( lib = 'openblas', fragment = blasCode, uselib_store = 'CBLAS' ):
			conf.define( 'HAVE_CBLAS', 1 ) ;
		if conf.check_cxx(
			lib = 'lapack',
			uselib_store = 'LAPACK'
		):
			conf.define( 'HAVE_LAPACK', 1 )
		if conf.check_cxx( header_name = 'sys/time.h', fragment = '#include "sys/time.h"\nint main(int argc, char** argv ) { struct timeval current_time ; gettimeofday( &current_time, 0 ) ; return 0 ; }' ):
			conf.define( 'HAVE_GETTIMEOFDAY', 1 )
			
def misc_configure( conf ) :
	conf.define( 'EIGEN_NO_DEBUG', 1 )

def get_cxxflags( variant_name ):
	cxxflags = [
		'-Wall',
		'-pedantic',
		'-Wno-long-long', # don't warn about the long long thing, it comes up in Eigen and Boost.
		#'-Wno-redeclared-class-member', # don't warn about class member redeclaration which comes up in Boost
		'-Wno-unused-local-typedefs' # warns in boost
	]
	if variant_name == 'default':
		cxxflags.extend( [ '-g' ] )
	elif variant_name == 'release':
		cxxflags.extend( [ '-O3' ])
	return cxxflags

def get_ldflags( variant_name ):
	import platform
	ldflags = []
	if variant_name == 'default' and platform.system() == 'Darwin':
		ldflags.extend( [ '-framework', 'CoreFoundation' ])
	if Options.options.static and platform.system() != 'Darwin':
		ldflags.extend( [ '-static', '-static-libgcc' ] )
	return ldflags

def create_variant( conf, variant_name ):
	variant = conf.env.copy()
	conf.set_env_name( variant_name, variant )
	variant.set_variant( variant_name )

def configure_variant( conf, variant_name, cxxflags = [], linkflags = [] ):
	cxxflags.extend( get_cxxflags( variant_name ))
	cxxflags.extend( [ '-I', variant_name ] )
	linkflags.extend( get_ldflags( variant_name ))
	
	conf.setenv( variant_name )
	conf.env[ 'CXXFLAGS' ] = cxxflags
	conf.env[ 'LINKFLAGS' ] = linkflags
	conf.write_config_header( 'config/config.hpp' )
#	conf.write_config_header( 'config.hpp' )
#	conf.write_config_header( 'genfile/config.hpp' )

#-----------------------------------
# BUILD
#-----------------------------------

def build( bld ):
	import Options

	bld(
                rule = """printf '#ifndef QCTOOL_REVISION_HPP\n#define QCTOOL_REVISION_HPP\nnamespace globals {\n\tchar const* qctool_version = \"%%s\" ;\n\tchar const* const qctool_revision = \"%%s\" ;\n}\n#endif\n' `echo """ + VERSION + "` `hg parents --template={node}` > ${TGT}""",
		#rule = """printf '#ifndef QCTOOL_VERSION_HPP\n#define QCTOOL_VERSION_HPP\nnamespace globals {\n\tchar const* const qctool_revision = \"%%s\" ;\n}\n#endif\n' `hg parents --template={node}` > ${TGT}""",
		always = True,
		target = "config/qctool_version_autogenerated.hpp",
		name = "qctool_version_autogenerated",
		uselib = "",
		#on_results = True
		on_results = False
	)
	
	for subdir in subdirs:	
		bld.add_subdirs( subdir )
	
	#---------------------
	# libs
	#---------------------

	bld.new_task_gen(
		features = 'cxx cstaticlib',
		target = 'gen-tools-lib',
		source = bld.glob( 'src/*.cpp' ),
		includes='./include ./genfile/include',
		uselib_local = 'string_utils statfile appcontext fputils worker genfile integration',
		uselib = 'BOOST BOOST_IOSTREAMS ZLIB BOOST_MATH BOOST_FILESYSTEM BOOST_SYSTEM MGL CBLAS CLAPACK MGL'
	)

	#---------------------
	# programs
	#---------------------
	USELIB = "BOOST BOOST_IOSTREAMS ZLIB BOOST_FILESYSTEM BOOST_SYSTEM BOOST_UNIT_TEST_FRAMEWORK BOOST_THREAD BOOST_TIMER BOOST_CHRONO PTHREAD CBLAS CLAPACK MGL RT"
	Components = 'SNPSummaryComponent SampleSummaryComponent HaplotypeFrequencyComponent RelatednessComponent SNPOutputComponent'
	create_app( bld, name='qctool', uselib = USELIB, uselib_local = Components + ' gen-tools-lib qcdb appcontext statfile appcontext string_utils fputils worker genfile metro qctool_version_autogenerated' )
	create_app( bld, name='hptest', uselib = USELIB, uselib_local = Components + ' gen-tools-lib qcdb appcontext statfile appcontext string_utils fputils worker genfile metro qctool_version_autogenerated' )
	create_app( bld, name='multifinemap', uselib = USELIB, uselib_local = Components + ' gen-tools-lib qcdb appcontext statfile appcontext string_utils fputils worker genfile metro qctool_version_autogenerated' )
	create_app( bld, name='inthinnerator', uselib = USELIB, uselib_local = Components + ' qctool_version_autogenerated gen-tools-lib genfile qcdb statfile appcontext' )

	if Options.options.all_targets:
		create_app( bld, name='overrep', uselib = 'BOOST_REGEX ' + USELIB, uselib_local = 'qctool_version_autogenerated qcdb appcontext gen-tools-lib genfile statfile' )
		create_app( bld, name='gen-grep', uselib = USELIB, uselib_local = 'gen-tools-lib genfile appcontext string_utils' )
		create_app( bld, name='inflation', uselib = USELIB, uselib_local = 'qctool_version_autogenerated appcontext genfile' )
		create_app( bld, name='selfmap', uselib = USELIB, uselib_local = 'SNPSummaryComponent qctool_version_autogenerated qcdb appcontext genfile' )
		create_app( bld, name='checkvcf', uselib = USELIB, uselib_local = 'SNPSummaryComponent qctool_version_autogenerated qcdb appcontext statfile genfile' )
		create_app( bld, name='binit', uselib = USELIB, uselib_local = 'qctool_version_autogenerated appcontext genfile' )
	
	# Build release variants as well.
	for obj in [] + bld.all_task_gen:
	    install_path = obj.install_path
	    new_obj = obj.clone('release')
	    new_obj.install_path = install_path
	    obj.install_path = None

	#---------------------
	# tests
	#---------------------
	# misc tests...

	# create_tests( bld, uselib = USELIB, components = Components )

def create_app( bld, name, uselib = '', uselib_local = '' ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target =  '%s_v%s' % ( name, VERSION ),
		source = [  'apps/' + name + '.cpp' ],
		includes='./ ./include ./genfile/include ./statfile/include',
		uselib_local = uselib_local,
		uselib = uselib + ' boost',
		install_path = os.path.join( bld.env[ 'PREFIX' ], 'bin' )
	)

def create_tests( bld, uselib = '', components = '' ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = 'test_qctool',
		source = bld.glob( 'test/*.cpp' ),
		uselib_local = 'gen-tools-lib genfile appcontext string_utils ' + components + ' boost-unit_test_framework',
		includes='./include ./genfile/include',
		uselib = uselib + ' BOOST_UNIT_TEST_FRAMEWORK',
		unit_test=1,
		install_path=None
	)

def create_benchmark( bld, name, uselib = '' ):
	bld.new_task_gen(
		features = 'cxx cprogram',
		target = name,
		source = [  'benchmarks/' + name + '.cpp' ],
		uselib_local = 'gen-tools-lib string_utils genfile',
		uselib = uselib + ' BOOST_UNIT_TEST_FRAMEWORK',
		includes='./include ./genfile/include',
		install_path=None
	)

#-----------------------------------
# CHECK
#-----------------------------------

def test(context):
	print( "Performing functional tests..." )
	import json
	working_dir = "release/test_data"
	qctool_executable = "build/release/qctool-%s" % VERSION
	json = json.loads( open( os.path.join( working_dir, "catalogue.json" ) ).read() )
	import sys
	sys.path.append( "release" )
	import Release.TestHarness
	harness = Release.TestHarness.TestHarness( qctool_executable, working_dir, json )
	harness.run()

	# UnitTest module no longer works in waf 1.5.18.  Need to update this.
	#import UnitTest
	#ut = UnitTest.unit_test()
	#ut.change_to_testfile_dir = True
	#ut.run()
	#ut.print_results()

def release( bld ):
	import sys
	import Options
	sys.path.append( "release" )
	import Release.TestHarness
	import Release.ReleaseBuilder
	
	if Options.options.all_targets:
		inthinnerator_executable = "build/release/inthinnerator_v%s" % VERSION
		builder = Release.ReleaseBuilder.ReleaseBuilder( "inthinnerator", VERSION, inthinnerator_executable )
		release = builder.build()
		print( "++ inthinnerator release tarball created in", release[ "release_tarball" ] )

	qctool_executable = "build/release/qctool_v%s" % VERSION
	builder = Release.ReleaseBuilder.ReleaseBuilder( APPNAME, VERSION, qctool_executable )
	release = builder.build()
	print( "++ qctool release tarball created in", release[ "release_tarball" ] )

	print( "++ Release tarball created in", release[ "release_tarball" ] )
