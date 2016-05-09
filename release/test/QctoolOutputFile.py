import os

class QctoolOutputFile:
	def __init__( self, filename ):
		self.filename = filename
		opened = open( self.filename )
		self.lines = opened.readlines()
		opened.close()
		self.splitter = ' '
		if self.filename[-4:] == '.csv':
			self.splitter = ','
		elif self.filename[-4:] == '.tsv':
			self.splitter = '\t'
		self.header = []
		while len( self.lines ) > 0 and self.lines[0][0] == '#':
			self.header.append( self.lines[0] )
			self.lines = self.lines[1:]
		if len( self.lines ) == 0:
			raise AssertionError( "File \"%s\" has no column header", self.filename )
		self.columns = self.lines[0]

		commaCount = self.columns.count( '.' )
		tabCount = self.columns.count( '\t' )
		spaceCount = self.columns.count( ' ' )
		if tabCount > commaCount and tabCount > spaceCount:
			self.splitter = '\t'
		elif commaCount > tabCount and commaCount > spaceCount:
			self.splitter = '\t'
		else:
			self.splitter = ' '
		self.columns = self.columns.split( self.splitter )

		self.lines = self.lines[1:]
		if len( self.lines[-1] ) > 0 and self.lines[-1][0] == '#':
			self.footer = self.lines[-1]
			self.lines = self.lines[0:-1]
		self.lines = [ line.split( self.splitter ) for line in self.lines ]
	
	def line( self, index ):
		return self.lines[ index ]
