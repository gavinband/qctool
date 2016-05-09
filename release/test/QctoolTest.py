import os
import re
from TextFile import TextFile
from QctoolOutputFile import QctoolOutputFile
from robot.api import logger

class QctoolTest:
	ROBOT_LIBRARY_SCOPE = 'TEST CASE'

	def has_valid_metadata( self, filename ):
		file = open( filename )
		lines = file.readlines( 1000 ) ;
		file.close()
		count = 0
		if lines[0][0:11] != '# Analysis:':
			raise AssertionError( 'File "%s" has invalid Analysis line ("%s")' % ( filename, lines[0].strip() ) )
		if lines[1][0:11] != '#  started:':
			raise AssertionError( 'File "%s" has invalid started line ("%s")' % ( filename, lines[1][0:11] ) )
		if lines[3][0:22] != '# Analysis properties:':
			raise AssertionError( 'File "%s" has invalid Analysis properties line("%s")' % ( filename, lines[3].strip() ) )

	def has_column_header( self, filename ):
		"""Checks that output file has header lines including a column name line"""
		file = open( filename )
		lines = file.readlines( 1000 ) ;
		file.close()
		count = 0
		while len( lines ) > 0 and len(lines[0]) > 0 and lines[0][0] == '#':
			count += 1
			lines = lines[1:]

		if len( lines ) == 0:
			raise AssertionError( "File \"%s\" has no column header line" % filename )
		elts = lines[0].split( '\t' )
	
		if len(elts) < 2:
			logger.info( '-'.join( elts ))
			raise AssertionError( "File \"%s\" has malformed column header" % filename )
	
	def has_footer( self, filename ):
		file = open( filename )
		lines = file.readlines() ;
		file.close()
		expected = '# Completed successfully at'
		if lines[-1][0:len( expected )] != expected:
			raise AssertionError( 'File ("%s") has invalid footer' % filename )

	def has_column( self, filename, column ):
		data = QctoolOutputFile( filename )
		if not column in data.columns:
			raise AssertionError( "File \"%s\" does not contain expected column \"%s\"" % ( filename, column ) )

	def has_column_matching( self, filename, columnRegexp ):
		columnRe = re.compile( columnRegexp )
		data = QctoolOutputFile( filename )
		for column in data.columns:
			if re.match( columnRe, column ):
				return
		raise AssertionError( "No column in file \"%s\" matches the regular expression \"%s\": %s." % ( filename, columnRegexp, ' '.join( data.columns ) ) )
	
	def matches_linecount( self, filename, other_filename, other_skip ):
		data = QctoolOutputFile( filename )
		other_file = open( other_filename )
		n = len( other_file.readlines() ) - int( other_skip )
		other_file.close()
		if len( data.lines ) != n:
			raise AssertionError( "File \"%s\" has wrong number of lines (%d instead of %d)" % ( filename, len( data.lines ), n ) )
	
	def columns_match( self, referenceFileName, referenceFileColumn, targetFileName, targetFileColumn, decimalPlaces = None ):
		referenceFileColumnRegexp = re.compile( referenceFileColumn )
		targetFileColumnRegexp = re.compile( targetFileColumn )
		referenceFile = QctoolOutputFile( referenceFileName )
		targetFile = QctoolOutputFile( targetFileName )
		referenceFileIndex = [ i for i, item in enumerate( referenceFile.columns ) if re.match( referenceFileColumnRegexp, item ) ]
		targetFileIndex = [ i for i, item in enumerate( targetFile.columns ) if re.match( targetFileColumnRegexp, item ) ]
		if len( referenceFileIndex ) == 0:
			raise AssertionError( 'File "%s" has no column matching "%s"' % ( referenceFileName, referenceFileColumn ))
		if len( targetFileIndex ) == 0:
			raise AssertionError( 'File "%s" has no column matching "%s"' % ( targetFileName, targetFileColumn ))

		referenceFileIndex = referenceFileIndex[0]
		targetFileIndex = targetFileIndex[0]

		floatFormat = None
		if decimalPlaces is not None:
			if len( decimalPlaces ) == 3 and decimalPlaces[1:] == 'sd':
				floatFormat = "%%.%sg" % decimalPlaces[0]
			else:
				floatFormat = "%%.%sf" % decimalPlaces

		for i in range( 0, len( referenceFile.lines ) ):
			referenceElt = referenceFile.lines[i][ referenceFileIndex ]
			targetElt = targetFile.lines[i][ targetFileIndex ]
			if not self.values_match( referenceElt, targetElt, floatFormat, allowNewNonMissingValues = True ):
				raise AssertionError(
					'Column "%s" of file "%s" does not match column "%s" of file "%s" in the %dth entry ("%s", "%s")'
					% ( referenceFileColumn, referenceFileName, targetFileColumn, targetFileName, i, referenceFile.lines[i][ referenceFileIndex ], targetFile.lines[i][ targetFileIndex ] )
				)

	# Test whether two values match up to some choices of behaviour
	# Values equal to 'NA' are treated as missing values.  If allowNewNonMissingValues=True
	# then the targetElt may be non-missing even if referenceElt is missing.  Otherwise missingness
	# must be the same for both values.
	# If both values are non-missing, if decimalPlaces is not None then values are treated as numbers
	# and compared to be equal to the given number of decimal places.  Otherwise they are treated
	# as strings and compared literally.
	def values_match( self, referenceElt, targetElt, floatFormat = None, allowNewNonMissingValues = False ):
		return (
			( referenceElt == 'NA' and targetElt == 'NA' )
			or
			( allowNewNonMissingValues and referenceElt == 'NA' and targetElt != 'NA' )
			or
			( floatFormat is None and referenceElt == targetElt )
			or
			( floatFormat is not None and referenceElt != 'NA' and targetElt != 'NA' and (floatFormat % float(referenceElt)) == (floatFormat % float(targetElt)) )
		)
	
