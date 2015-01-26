import os
import subprocess
import tempfile
import filecmp

class SingleTest:
	def __init__( self, command, output_dir, result ):
		self.command = command
		self.result = result
		self.output_dir = output_dir
		self.directory_comparison = None
		self.file_comparison = None
		if not os.path.exists( output_dir ) or not os.path.isdir( output_dir ):
			raise Exception( "SingleTest: The output dir \"" + output_dir + "\" does not exist." )
		if ( not os.path.exists( result ) ):
			raise Exception( "SingleTest: The specified result \"" + result + "\" does not exist." )

	def run( self ):
		subprocess.check_call( self.command.split( " " ), stderr = open( "/dev/null", 'w' ) )
		expected_result_dir = self.extract_expected_result()
		self.compare( expected_result_dir, self.output_dir )
	
	def extract_expected_result( self ):
		if os.path.isfile( self.result ) and ( self.result.endswith( ".tgz" ) or self.result.endswith( ".tar.gz" ) ):
			# Extract the expected result tarball
			tmpdir = tempfile.mkdtemp( prefix="qctool-expected" )
			subprocess.check_call( [ "tar", "-C", tmpdir, "-xzf", self.result ] )
			return tmpdir
		else:
			return self.result

	def compare( self, expected_dir, output_dir ):
		print "...and comparing to", expected_dir

		# Perform a comparison
		F = filecmp.dircmp( expected_dir, output_dir )
		self.directory_comparison = F
		common_files = F.common_files
		F = filecmp.cmpfiles( expected_dir, output_dir, common_files, shallow = False )
		self.file_comparison = { "match": F[0], "mismatch": F[1], "errors": F[2] }

		if len( self.directory_comparison.left_only ) > 0 or len( self.directory_comparison.right_only ) > 0 or len( self.directory_comparison.common_funny ) > 0:
			print '!! For command', self.command, 'expected result and actual result differ:'
			self.directory_comparison.report()
			raise Exception( "Expected result differs from actual result" )

		if ( len( self.file_comparison[ "mismatch" ] ) + len( self.file_comparison[ "errors" ] ) ) != 0:
			print '!! For command', self.command, 'expected result and actual result differ:'
			print "Mismatches:", " ".join( self.file_comparison[ "mismatch" ] ), "Errors:", " ".join( self.file_comparison[ "errors" ] )
			raise Exception( "Expected result differs from actual result" )

	def get_comparison( self ):
		return { "directory_comparison": self.directory_comparison, "file_comparison": self.file_comparison }
