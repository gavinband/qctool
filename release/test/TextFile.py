import os

class TextFile:
	def __init__( self, filename ):
		self.filename = filename
		opened = open( self.filename )
		self.lines = opened.readlines()
		opened.close()
		self.header = []
		while len( self.lines ) > 0 and self.lines[0][0] == '#':
			self.header.append( self.lines[0] )
			self.lines = self.lines[1:]
		if len( self.lines ) == 0:
			raise AssertionError( "File \"%s\" has no column header", self.filename )
		self.columns = self.lines[0].split( ' ' )
		self.lines = self.lines[1:]
		if len( self.lines[-1] ) > 0 and self.lines[-1][0] == '#':
			self.footer = self.lines[-1]
			self.lines = self.lines[0:-1]
		self.lines = [ line.split( ' ') for line in self.lines ]
	
	def line( self, index ):
		return self.lines[ index ]
