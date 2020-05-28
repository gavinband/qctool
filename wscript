import os.path
import glob

srcdir="."
APPNAME = "qctool"
VERSION = "2.1.2"

variants = {
	"release": {
		"command": "{command}",
		"cxxflags": [ '-O3' ],
		"linkflags": []
	},
	"debug": {
		"command": "{command}_debug",
		"cxxflags": [ '-g' ],
		"linkflags": []
	}
} ;

subdirs = [
	'genfile', 'statfile', 'appcontext',
	'worker', 
	'3rd_party', 'components', 'qcdb', 'metro',
	'apps'
]

def options( opt ):
	opt.load( 'compiler_cxx' )
	opt.load( 'compiler_c' )
	opt.add_option( "--all_targets", action='store_true', default=False, help='Create all targets, not just qctool.')
	opt.add_option( "--vectorise", action='store_true', default=False, help='Use aggresive vectorisation options for maximum performance.')
	opt.add_option( '--variant', action='store', default='release', help='set the variant name' )
	
#-----------------------------------
# CONFIGURE
#-----------------------------------

def configure( cfg ):
	configure_variant( cfg, 'release' )
	configure_variant( cfg, 'debug' )

def configure_variant( cfg, variant ):
	cfg.setenv( variant )
	cfg.env[ 'VERSION' ] = VERSION ;
	
	cfg.load( 'compiler_c' )
	cfg.load( 'compiler_cxx' )

	print( "Using prefix\t\t\t\t :", cfg.env[ 'PREFIX' ] )

	cxxflags = cfg.env[ 'CXXFLAGS' ] 
	linkflags = cfg.env[ 'LINKFLAGS' ]

	# We now require C++11.
	if cfg.check( cxxflags = '-std=c++11' ):
		cxxflags.append( '-std=c++11' )
		cxxflags.append( '-Wno-deprecated-declarations' )
	else:
		raise( "QCTOOL now requires a C++11-enabled compiler." )

	cxxflags.extend(
		[
			"-mavx", "-mssse3", "-msse2", "-msse4.1", "-msse4.2",
			'-Wall',
			'-pedantic',
			'-Wno-long-long', # suppress warnings in Eigen and Boost.
			# The following are to suppress many warnings in boost
			'-Wno-unused-local-typedefs',
			'-Wno-c++11-long-long',
			'-Wno-keyword-macro',
			'-Wno-unused-const-variable',
			'-Wno-deprecated-register',
			'-Wno-unused-function',
			'-Wno-redeclared-class-member'
		]
	)

	if cfg.options.vectorise:
		# These are disabled by default as not always safe.
		for flag in [ '-msse2', '-mavx', '-mssse3']:
			if cfg.check( cxxflags = flag ):
				cxxflags.append( flag )
	
	cxxflags.extend( variants[ variant ][ "cxxflags" ] )
	linkflags.extend( variants[ variant ][ "linkflags" ] )
	
	configure_blas( cfg )
	configure_time( cfg )
	import platform
	if platform.system() == 'Darwin':
		configure_darwin( cfg, cxxflags, linkflags )
	else:
		configure_linux( cfg, cxxflags, linkflags )

	# Now check for libraries
	if check_cxx( cfg, lib = 'm', uselib_store = 'M' ):
		cfg.define( 'HAVE_M', 1 )
	if check_cxx( cfg, lib='z', uselib_store='ZLIB' ):
		cfg.define( 'HAVE_ZLIB', 1 )
	# disable support for some things
	cfg.define( 'HAVE_BZIP2', 0 )
	cfg.define( 'HAVE_MGL', 0 )
	cfg.define( 'HAVE_CAIRO', 0 )
	# sqlite3, eigen, zstd are part of this repo
	cfg.define( 'HAVE_SQLITE3', 1 )
	cfg.define( 'HAVE_EIGEN', 1 )
	cfg.define( 'HAVE_ZSTD', 1 )
	# Boost libs are now contained in the repo
	cfg.define( "HAVE_BOOST_IOSTREAMS", 1 )
	cfg.define( "HAVE_BOOST_FILESYSTEM", 1 )
	cfg.define( "HAVE_BOOST_SYSTEM", 1 )
	cfg.define( "HAVE_BOOST_THREAD", 1 )
	cfg.define( "HAVE_BOOST_DATE_TIME", 1 )
	cfg.define( "HAVE_BOOST_UNIT_TEST_FRAMEWORK", 1 )
	cfg.define( "HAVE_BOOST_TIMER", 1 )
	cfg.define( "HAVE_BOOST_REGEX", 1 )
	cfg.define( "HAVE_BOOST_MATH", 1 )
	cfg.define( "HAVE_BOOST_FUNCTION", 1 )
	cfg.define( "HAVE_BOOST_SPIRIT", 1 )
	cfg.define( "BOOST_FILESYSTEM_VERSION", 3 )
	cfg.define( 'EIGEN_NO_DEBUG', 1 )

	cxxflags.append( "-I." )
	cfg.env[ 'CXXFLAGS' ] = cxxflags
	cfg.env[ 'LINKFLAGS' ] = linkflags
	cfg.write_config_header( '%s/config/config.hpp' % variant )

def check_cxx( cfg, **kwargs ):
	try:
		cfg.check_cxx( **kwargs )
		return True
	except:
		return False
		

def configure_blas( cfg ):
	import platform
	if platform.system() == 'Darwin':
		if cfg.check_cxx(
			lib = 'cblas',
			fragment = '#include "cblas.h"\nint main() {}',
			cxxflags = '-I/System/Library/Frameworks/vecLib.framework/Headers',
			uselib_store = 'CBLAS'
		):
			cfg.define( 'HAVE_CBLAS', 1 )
		if check_cxx(
			cfg,
			header_name = 'clapack.h',
			lib = 'clapack',
			cxxflags = '-I/System/Library/Frameworks/vecLib.framework/Headers',
			uselib_store = 'CLAPACK'
		):
			cfg.define( 'HAVE_CLAPACK', 1 )
			cfg.define( 'HAVE_LAPACK', 1 )
	else:
		blasCode = """#include "cblas.h"
int main() {
	float v[10] ;
	float result = cblas_snrm2( 1, v, 1 ) ;
	return int(result) ;
}
"""
		if check_cxx( cfg, lib = 'cblas', fragment = blasCode, uselib_store = 'CBLAS' ):
			cfg.define( 'HAVE_CBLAS', 1 ) ;
		elif check_cxx( cfg, lib = 'blas', fragment = blasCode, uselib_store = 'CBLAS' ):
			cfg.define( 'HAVE_CBLAS', 1 ) ;
		elif check_cxx( cfg, lib = 'openblas', fragment = blasCode, uselib_store = 'CBLAS' ):
			cfg.define( 'HAVE_CBLAS', 1 ) ;
		if check_cxx(
			cfg,
			lib = 'lapack',
			uselib_store = 'LAPACK'
		):
			cfg.define( 'HAVE_LAPACK', 1 )

def configure_time( cfg ):
	if check_cxx( cfg, header_name='mach/mach_time.h', uselib_store = 'MACH_TIME' ):
		cfg.define( 'HAVE_MACH_TIME', 1 )
	if check_cxx( cfg, header_name = 'sys/time.h', fragment = '#include "sys/time.h"\nint main(int argc, char** argv ) { struct timeval current_time ; gettimeofday( &current_time, 0 ) ; return 0 ; }' ):
		cfg.define( 'HAVE_GETTIMEOFDAY', 1 )

def configure_darwin( cfg, cxxflags, linkflags ):
	linkflags.extend( [ '-framework', 'CoreFoundation' ])

def configure_linux( cfg, cxxflags, linkflags ):
    check_cxx( cfg, lib='rt', uselib_store='rt', msg = 'rt' )
    check_cxx( cfg, lib='pthread', uselib_store='pthread', msg = 'pthread' )
    check_cxx( cfg, lib='dl', uselib_store='dl', msg = 'dl' )

# This piece of code comes from
# https://gitlab.com/ita1024/waf/-/blob/master/demos/variants/wscript
# It makes the --variant command-line option work.	
def init(ctx):
    from waflib.Options import options
    from waflib.Build import BuildContext, CleanContext, InstallContext, UninstallContext
    from waflib.Configure import ConfigurationContext
    for y in (BuildContext, CleanContext, InstallContext, UninstallContext, ConfigurationContext):
        name = y.__name__.replace('Context','').lower()
        class tmp(y):
            cmd = name
            variant = options.variant

#-----------------------------------
# BUILD
#-----------------------------------

def create_app( bld, name, use = '' ):
	bld.program(
		target =  '%s_v%s' % ( name, VERSION ),
		source = 'apps/' + name + '.cpp',
		includes = './ include',
		use = use,
		install_path = os.path.join( bld.env[ 'PREFIX' ], 'bin' )
	)

def build( bld ):
	print( "Building %s variant..." % bld.variant )
	bld.env[ 'VERSION' ] = VERSION
	autogenerate_version( bld )

	bld.stlib(
		target = "qctool-lib",
		source = bld.path.ant_glob( "src/*.cpp" ),
		includes = "./include",
		export_includes = "./include",
		use = "genfile statfile appcontext"
	)
	
	for subdir in subdirs:	
		bld.recurse( subdir )
	
	
def autogenerate_version( bld ):
	bld(
		target = "config/version_autogenerated.hpp",
		rule = """printf '#ifndef QCTOOL_REVISION_HPP
#define QCTOOL_REVISION_HPP
namespace globals {
  char const* qctool_version = \"%%s\" ;
  char const* const qctool_revision = \"%%s\" ;
}
#endif
' \
`echo """ + VERSION + "` \
`fossil status | grep checkout | sed -e 's/checkout: *//' | cut -d' ' -f1` \
> ${TGT}""",
		name = "version_autogenerated",
		ext_out = [ '.hpp' ],
		always = True
	)

