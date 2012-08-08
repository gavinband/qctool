import os, os.path
import tempfile
import hashlib
import subprocess
import pprint
from . import SingleTest

class TestHarness:
	def __init__( self, qctool_exe, working_dir, data ):
		self.working_dir = working_dir
		self.qctool_exe = qctool_exe
		self.setup_tarball = os.path.join( working_dir, data[ "setup" ] )
		self.tests = data[ "tests" ]
		for i in range( 0, len( self.tests ) ):
			if not "result_tarball" in self.tests[i]:
				self.tests[i][ "result_tarball" ] = os.path.join( self.working_dir, hashlib.md5( self.tests[i][ "command" ] ).hexdigest() + ".tgz" )
			if not "name" in self.tests[i]:
				self.tests[i][ "name" ] = self.tests[i][ "command" ]
		print "Tests are:"
		pprint.pprint( [ test[ "name" ] for test in self.tests ] )
		
		if not os.path.exists( qctool_exe ) and os.path.isfile( qctool_exe ):
			raise Exception( "The qctool executable \"" + qctool_exe + "\" must exist and be an executable." )
		if not os.path.exists( working_dir ) and os.path.isdir( working_dir ):
			raise Exception( "The working dir \"" + working_dir + "\" must exist and be a directory." )
		self.m_success = None
		
	def run( self ):
		working_dir = self.setup()
		
		test_results = { "succeeded": [], "failed": [] }

		for i in range( 0, len( self.tests )):
			command = self.tests[i][ "command" ]
			output_dir = tempfile.mkdtemp( prefix="qctool-test-" )
			command = command.replace( '%INPUT_DIR%', working_dir )
			command = command.replace( '%OUTPUT_DIR%', output_dir )
			print "Running command: \"%s\"" % ( os.path.basename( self.qctool_exe ) + " " + command )
			
			test = SingleTest.SingleTest( self.qctool_exe + " " + command, output_dir, os.path.join( self.working_dir, self.tests[i][ "result_tarball" ] ))
			try:
				test.run()
				test_results[ "succeeded" ].append( self.tests[i][ "name" ] )
			except Exception, e:
				print "!! Testing of command \"" + self.tests[i][ "name" ] + "\" failed!"
				print e
				test_results[ "failed" ].append( { "name": self.tests[i][ "name" ], "comparison": test.get_comparison(), "command": command } )
		
		print "---- test summary ----"
		print "succeeded:     ", len( test_results[ "succeeded" ] )
		print "failed:        ", len( test_results[ "failed" ] )
			
		if len( test_results[ "failed" ] ) > 0:
			print "Tests failed for the following commands:"
			for i in range( 0, len( test_results[ "failed" ] )):
				pprint.pprint( test_results[ "failed" ][i] )
				print "   ( expected", self.tests[i][ "result_tarball" ], ")"
			self.m_success = False
		else:
			self.m_success = True
			
	def setup( self ):
		working_dir = tempfile.mkdtemp( prefix="qctool-test-" )
		print "Extracting", self.setup_tarball, "to", working_dir
		subprocess.check_call( [ "tar", "-C", working_dir, "-xzf", self.setup_tarball ] )
		return working_dir

	def success( self ):
		return self.m_success