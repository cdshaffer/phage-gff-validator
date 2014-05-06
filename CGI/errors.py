#!/usr/bin/env python

# Author: Paul Lee <pfleewustl@gmail.com>

"""
A file validator for the Phage Hunters class taught at Washington University in
St. Louis. Intended for use with gff_validator.html, SaveFile.cgi, and post_read.html.

Exception classes:

- 'Validation Error': basic exception class. Used to catch unknown errors. 
- 'Format Error': catches errors dealing with line format.
- 'Line Error': catches errors dealing with individuel components of a line.
- 'Biology Error': catches errors dealing with the described genes.
"""

class ValidationError(Exception):
	def __init__(self, code, message = ""):
        
		_dict = dict()
		_dict['0000'] = "Usage Error: code not used. Use printDict() to see dictionary of errors"
		_dict['1000'] = "Validation Error: unknown"
        
		self.code = str(code)
		self._dict = _dict
		self.message = message
        
		super(ValidationError, self).__init__(code, message)
        
	def printDict(self):
		for k, v in self._dict.iteritems():
			print k + ":  " + v
            
	def returnError(self):
		try:
			return self._dict[self.code]
		except:
			return self._dict['0000']


class FormatError(ValidationError):
	def __init__(self, code, message = ""):
        
		super(FormatError, self).__init__(code, message)
        
		self._dict['0100'] = "Format Error: unknown."
		self._dict['0200'] = "Format Error: lines merged."
		self._dict['0300'] = "Format Error: tab deliminated."
		self._dict['0400'] = "Format Error: multiple contigs. There should only be one contig line per file."
		self._dict['0500'] = "Format Error: each set of coordinates must have a line for each type. Types are seen at the third component. Default types are [gene, mRNA, exon]"
		self._dict['0600'] = "Format Error: ensure only one sequence in fasta file."
        

class LineError(ValidationError):
	def __init__(self, code, message = ""):
        
		super(LineError, self).__init__(code, message)

		self._dict['0001'] = "Line Error: unknown"
		self._dict['0002'] = "Line Error: 1st component = restrict characters used to a-Z/0-9/./=/;"
		self._dict['0003'] = "Line Error: 2nd component = restrict characters used to a-Z/0-9/./=/;"
		self._dict['0004'] = "Line Error: 3rd component = type not found in types. Default is [gene, mRNA, exon]."
		self._dict['0005'] = "Line Error: 4th component = 1st coordinate must be positive."
		self._dict['0006'] = "Line Error: 5th component = 2nd coordinate must be greater than the 1st coordinate."
		self._dict['0007'] = "Line Error: 6th component = must be a number or '.'"
		self._dict['0008'] = "Line Error: 7th component = must be '+' or '-'"
		self._dict['0009'] = "Line Error: 8th component = must be '.'"
		self._dict['0011'] = "Line Error: 9th component = restrict characters used to a-Z/0-9/./=/;"
		self._dict['0012'] = "Line Error: 9th component = ID already used. Each ID must be unique."
		self._dict['0013'] = "Line Error: 9th component = must have an ID and a name"
		self._dict['0021'] = "Line Error: 9th component = bad sort or a gene does not have a line for each type. Default types are [gene, mRNA, exon]."
		self._dict['0022'] = "Line Error: 9th component = restrict characters used to a-Z/0-9/./=/;"
		self._dict['0023'] = "Line Error: 9th component = ID already used. Each ID must be unique."
		self._dict['0024'] = "Line Error: 9th component = last type, must at least have a Parent."
		self._dict['0025'] = "Line Error: 9th component = must have an ID and a Parent."
        
        
class BiologyError(ValidationError):
	def __init__(self, code, message = ""):
        
		super(BiologyError, self).__init__(code, message)
        
		self._dict['0010'] = "Biology Error: unknown"
		self._dict['0020'] = "Biology Error: Incorrect start codon. Must start with M"
		self._dict['0030'] = "Biology Error: No stop codon detected."
		self._dict['0040'] = "Biology Error: Internal stop codons"
		self._dict['0050'] = "Biology Error: total nucleotide count not divisible by three."    