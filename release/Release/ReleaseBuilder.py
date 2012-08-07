class ReleaseBuilder:
	def __init__( self, APPNAME, VERSION, executable ):
		self.APPNAME = APPNAME
		self.executable = executable
		self.VERSION = VERSION

	def build( self ):
		import os, tempfile, shutil, subprocess
		tempdir = tempfile.mkdtemp()
		import platform
		if platform.system() == 'Darwin':
			release_stub = '%s_v%s-osx' % ( self.APPNAME, self.VERSION )
		elif platform.system() == 'Linux':
			release_stub = '%s_v%s-linux-%s' % ( self.APPNAME, self.VERSION, platform.machine() )
		release_dir = os.path.join( tempdir, release_stub )
		os.mkdir( release_dir )
		
		target_executable = "%s/%s" % ( release_dir, self.APPNAME )
		
		shutil.copyfile( self.executable, target_executable )
		shutil.copymode( self.executable, target_executable )
		shutil.copyfile( "LICENSE_1_0.txt", "%s/LICENSE_1_0.txt" % release_dir )
		shutil.copyfile( "CHANGELOG.txt", "%s/CHANGELOG.txt" % release_dir )
		target_tarball =  '%s/%s.tgz' % ( tempdir, release_stub )
		process = subprocess.Popen( [ 'tar', '-czf', target_tarball, release_stub ], cwd = tempdir )
		process.wait()
		print 'Created %s release tarball in "%s"' % ( self.APPNAME, target_tarball )
		print "Contents are:"
		print subprocess.Popen( [ 'tar', '-tzf', target_tarball ], stdout = subprocess.PIPE ).communicate()[0]
		return { "release_executable": target_executable, "release_tarball": target_tarball }
