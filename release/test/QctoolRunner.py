import os, sys
import subprocess
import tempfile
import json
from robot.api import logger
from StringIO import StringIO

class QctoolRunner:
	def __init__( self, qctool ):
		self.qctool = qctool
		self.command = 'snp-stats'
		self.gdata = [ 'data/test.gen' ]
		self.sdata = [ 'data/test.sample' ]
		self.log = None
		self.outputFilename = None

	def set_data( self, dataSpec ):
		spec = json.load( StringIO( dataSpec ))
		self.gdata = spec.get( "g" )
		self.sdata = spec.get( "s" )
	
	def run_qctool( self, command, expected, qctool = None ):
		if qctool is None:
			qctool = self.qctool
		filenames = {}
		temp = open( 'temp', 'a' )
		print >> temp, command
		values = json.load( StringIO( command ))
		gdata = values.get( 'g', self.gdata )
		sdata = values.get( 's', self.sdata )
		ogdata = values.get( 'og', None )
		if ogdata == '<tmp>':
			filenames[ 'og' ] = tempfile.mkstemp()[1]
			ogdata = filenames[ 'og' ]
		osdata = values.get( 'os', None )
		if osdata == '<tmp>':
			filenames[ 'os' ] = tempfile.mkstemp()[1]
			osdata = filenames[ 'os' ]
		cmd = [ qctool ]
		logger.debug( gdata, sdata )
		for i in range( 0, len(gdata) ):
			cmd += [ "-g", gdata[i] ]
		for i in range( 0, len(sdata) ):
			cmd += [ "-s", sdata[i] ]
		
		for option in [ 'vcf-genotype-field', 'incl-samples', 'excl-samples', 'osnp', 'osample' ]:
			if values.get( option, None ) is not None:
				value = values[ option ]
				if value == '<tmp>':
					filenames[ option ] = tempfile.mkstemp()[1]
					value = filenames[ option ]
				cmd.extend( [ '-%s' % option, value ] )

		for option in [ 'snp-stats', 'sample-stats' ]:
			if values.get( option, None ) is not None:
				cmd.extend( [ '-%s' % option ] )
		
		if ogdata is not None:
			cmd.extend( [ '-og', ogdata ])
		if osdata is not None:
			cmd.extend( [ '-os', osdata ])
		
		print >> temp, ' '.join( cmd )
		try:
			subprocess.check_call( cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE )
			if expected != 'ok':
				raise Exception( "Expected failure: %s" % command )
		except subprocess.CalledProcessError, e:
			if expected == 'ok':
				raise Exception( "Expected success: %s" % command )
		finally:
			temp.close()
		logger.info( "Complete run" )
		return filenames

	def output_filename( self ):
		return self.outputFilename
