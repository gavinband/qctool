#!/usr/bin/env python3

import glob
import markdown
import os
import argparse

parser = argparse.ArgumentParser( description = """Build a documentation website by processing .content.html files into .html files.
A header and footer is added to each file, and the contents are processed through python-markdown.""" )
parser.add_argument(
	'--path', '-p',
	help = 'Path to operate on'
)
parser.add_argument(
	'--header', 
	help = 'File containing header text, relative to --path',
	default = 'header.html'
)
parser.add_argument(
	'--footer',
	help = 'File containing footer text, relative to --path',
	default = 'footer.html'
)

args = parser.parse_args()

print( "Building docs in %s..." % args.path )

contentPath = os.path.join( args.path, 'content' )

headerFilename = "%s/header.html" % contentPath
footerFilename = "%s/footer.html" % contentPath
print( "++ Reading header from %s..." % headerFilename )
header = open( headerFilename, 'r' ).read()
print( "++ Reading footer from %s..." % headerFilename )
footer = open( footerFilename, 'r' ).read()

initialDepth = contentPath.count(os.sep)
total = 0
for directory, dnames, filenames in os.walk( contentPath ):
	directoryParts = directory.split( os.sep )
	outputDirectory = os.sep.join( [ directoryParts[0] ] + directoryParts[2:] )
	directoryLevel = (directory.count( os.sep ) - initialDepth)
	for filename in filenames:
		print( "++ Looking at %s..." % filename )
		if filename.endswith( ".content.html" ):
			total = total + 1
			print( "++ Reading content from %s..." % os.path.join( directory, filename ))
			content = open( os.path.join( directory, filename ), 'r' ).read() ;
			content = (
				header.replace( "[ROOT]/", "../" * directoryLevel ) + "\n"
				+ markdown.markdown(
					content.replace( "[ROOT]/", "../" * directoryLevel ),
					extensions=['markdown.extensions.md_in_html', 'mdx_math', 'markdown.extensions.fenced_code', 'tables' ],
					extension_configs = {
						"mdx_math": { "enable_dollar_delimiter": True }
					}
				) + "\n"
				+ footer.replace( "[ROOT]/", "../" * directoryLevel )
			)

			outputFilename = os.path.join( outputDirectory, filename.replace( ".content.html", ".html" ))
			print( "++ Writing %s..." % outputFilename )
			outputFile = open( outputFilename, 'w' )
			outputFile.write( content )
			outputFile.close()

print( "++ Success, formatted a total of %d files." % total )
print( "++ Thanks for using build.py" )
