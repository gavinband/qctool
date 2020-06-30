import os, sys
import subprocess
import tempfile
import json
from robot.api import logger
try:
	# Python 2
	from StringIO import StringIO
except ImportError:
	# Python 3
	from io import StringIO

class QctoolRunner:
	def __init__( self, qctool ):
		self.qctool = qctool
		self.command = 'snp-stats'
		self.gdata = [ 'data/test.gen' ]
		self.sdata = None
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
		print( command, file = temp )
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
		if sdata is not None:
			for i in range( 0, len(sdata) ):
				cmd += [ "-s", sdata[i] ]
		
		for option in [
			'vcf-genotype-field',
			'incl-positions', 'excl-positions',
			'incl-snpids', 'excl-snpids', 'incl-rsids', 'excl-rsids',
			'incl-samples', 'incl-samples-where', 'excl-samples', 'excl-samples-where',
			'osnp', 'osample',
			'ofiletype', 'filetype',
			'annotate-bed3', 'annotate-bed4',
			'merge-in', 'assume-chromosome',
			'threshold', 'reorder', 'infer-ploidy-from',
			'strand', 'flip-to-match-allele',
			'bgen-bits'
		]:
			if values.get( option, None ) is not None:
				cmd.extend( [ '-%s' % option ] )
				valueList = values[ option ]
				if not isinstance( valueList, list ):
					valueList = [ valueList ]
				for i in range( 0, len(valueList)):
					if valueList[i] == '<tmp>':
						filenames[ option ] = tempfile.mkstemp()[1]
						valueList[i] = filenames[ option ]
					cmd.append( valueList[i] )

		for option in [ 'snp-stats', 'sample-stats' ]:
			if values.get( option, None ) is not None:
				cmd.extend( [ '-%s' % option ] )

		for option in [ 'force' ]:
			if option in values:
				cmd.extend( [ '-%s' % option ] )
		
		if ogdata is not None:
			cmd.extend( [ '-og', ogdata ])
		if osdata is not None:
			cmd.extend( [ '-os', osdata ])
		
		print( ' '.join( cmd ), file = temp )
		try:
			#subprocess.check_call( cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE )
			subprocess.check_call( cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL )
			if expected != 'ok':
				raise Exception( "Expected failure: %s" % command )
			print( 'ok', file = temp )
		except subprocess.CalledProcessError as e:
			if expected == 'ok':
				raise Exception( "Expected success: %s" % command )
		finally:
			temp.close()
		logger.info( "Complete run" )
		return filenames

	def output_filename( self ):
		return self.outputFilename
