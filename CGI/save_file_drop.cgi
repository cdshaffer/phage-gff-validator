#!/usr/bin/python

import cgi, os, sys, re
import cgitb
import gff_validator
import datetime
import time
import tempfile

MAX_FILE_SIZE = 100000

def main():

	status = "good"	
	incLine = False
	typeArr = []
	
	form = cgi.FieldStorage()
	keyList = form.keys()
	for key in keyList:

		if key == 'gffFile':
			gffItem = form[key]
			
		elif key == 'seqFile':
			seqItem = form[key]
			
		elif key == 'optionalErr':
			incLine = True
			
		elif key[0:4] == 'type':
			typeArr.append((str(key)[4:5], str(form[key].value)))
	
	if not typeArr[0][1]:
		typeArr = ['gene','mRNA','exon']
	else:
		typeArr = sortArray(typeArr)
		
	
	if gffItem.filename:
	
		gffItem.file.seek(0,2)
		filesize = gffItem.file.tell()
		gffItem.file.seek(0)
		
		
		if filesize > MAX_FILE_SIZE:
			sys.exit('GFF file to big.')
	
		fnA = os.path.basename(gffItem.filename.replace("\\", "/" ))
		dest_dir = '/Library/WebServer/tmp/'
		open(dest_dir + fnA, 'w').write(gffItem.file.read())
		os.rename(dest_dir + fnA, dest_dir + 'gffITEM.gff')
	else:
		status = "bad"

	if seqItem.filename:
	
		seqItem.file.seek(0,2)
		filesize = seqItem.file.tell()
		seqItem.file.seek(0)
		
		
		if filesize > MAX_FILE_SIZE:
			sys.exit('FASTA file to big.')

		
		fnB = os.path.basename(seqItem.filename.replace("\\", "/" ))
		dest_dir = '/Library/WebServer/tmp/'
		open(dest_dir + fnB, 'w').write(seqItem.file.read())
		os.rename(dest_dir + fnB, dest_dir + 'seqITEM.fasta')
		
		#open('/tmp/' + fnB, 'w').write(seqItem.file.read())	
	else:
		status = "bad"
	
	if status == "bad":
		print('ERROR: Problem reading either the gff or fasta file')
		sys.exit(-1)
	
	suf = '_DT_' + str(datetime.datetime.fromtimestamp(time.time())) +'.txt'
	suf.replace(" ","")
	dir = '/Library/WebServer/trash/'
	
	newErrors = tempfile.NamedTemporaryFile(suffix=suf, prefix='Errors_', dir=dir, delete=False)
	newSorted = tempfile.NamedTemporaryFile(suffix=suf, prefix='Sorted_', dir=dir, delete=False)
			
	gff_validator.main('/Library/WebServer/tmp/gffITEM.gff','/Library/WebServer/tmp/seqITEM.fasta', newErrors, newSorted, incLine, typeArr)
		
	newErrors.close()
	newSorted.close()
	
	item_E = "http://localhost/" + newErrors.name[19:]
	item_S = "http://localhost/" + newSorted.name[19:]
		
	new_html = '''
	<!DOCTYPE html>
	<html>
	
	<head>
		<title>Validation Results</title>
		<style type="text/css"></style>
	</head>
	
	<body>
	
		<p><a href="{item_E}">Errors</a></p>
		
		<p><a href="{item_S}">Sorted</a></p>
		
	</body>
	</html>
	'''
	
	print(new_html.format(**locals()))
	
	   
def sortArray(inArray):
	length = len(inArray)
	i = 1
	while i < length:
		 a = inArray[i-1]
		 b = inArray[i]
		 if a[0] > b[0]:
		 	inArray[i-1] = b
		 	inArray[i] = a
		 	if i > 1:
		 		i = i - 1
		 else:
		 	i = i + 1
	return inArray
	

try:
    print("Content-type: text/html\n\n")   # say generating html
    main() 
except:
    cgi.print_exception()                 # catch and print errors