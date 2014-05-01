#!/usr/bin/env python

# Author: Paul Lee <pfleewustl@gmail.com>

"""
A file validator for the Phage Hunters class taught at Washington University in
St. Louis. Intended for use with gff_validator.html, SaveFile.cgi, and post_read.html.

Exception classes:

- 'Validation Error': basic exception class. Used to catch unknown errors. 
- 'Format Error': catches errors dealing with line format.
- 'Line Error': catches errors dealing with individuel components of a line.
- 'Biology Error': catchers errors dealing with the genes provided.

Functions:

- 'main()': runs the module and outputs the errors found.
- 'sortGff3()': sorts the lines in the document based on line type.
- 'fileCheck()': checks each line in for proper format.
- 'charCheck()': checks a string for specific characters.
- 'geneCheck()': checks that the sequence given is a gene.
- 'fastaRead()': reads in a .fasta file to a string.
- 'translate()': reads a nucleotide sequence and outputs the protein sequence. 

How to Use This Module
======================
(see the individual classes, methods, and attributes for details.)

1. Input a file in .gff format and a file in .fasta format to main().
	- optional input: include input lines in output file.
	- optional input: type hierarchy

2. Retrieve output file from given directory.
	- output file is a basic .txt file.
	- NOTE: files will be overwritten each time the module is run. 
"""

import sys
import re
import os

"""
The main method of the module.

Parameters:
- 'gff': a text file, expected to be in .gff format.
- 'seq': a text file, expected to be in .fasta format.
- 'incLine': boolean indicating whether lines from the gff file should be included in the
				error output file.
- 'typeHier': a list indicating the types of each gff line and the order the types should
				be sorted in. The first type in the list will be ordered before the second, 
				the second type before the third, etc.
"""
def main(gff, seq, newErrors, newSorted, incLine=False, typeHier=['gene','mRNA','exon']):
		
	gff3_File = gff
	seq_File = seq
	seq = fastaRead(seq_File)
	
	global Errors
	Errors = []
    
	sorted_File = sortGff3(gff3_File, typeHier)    
	report = fileCheck(sorted_File[0], sorted_File[1], seq, typeHier)
    	
	outFile(sorted_File[0], sorted_File[1], newSorted, Errors, newErrors, incLine)
	return

"""
Sorts each line in the gff file according to the type hierarchy. The given type (ie 
'gene', 'mRNA', 'exon' must be identical to the types found in the file. 
Gene != gene, mrna != mRNA, etc.

Also checks lines for correct format:
- each line is tab delimited and has 9 components.
- each line has a type at the 3rd component.
- example line:
1stComp\t2ndComp\t3rdComp\t4thComp\t...

Parameters:
-'gff3_File': the uploaded gff file.
-'types': a list indicating the types of each gff line and the order the types should
				be sorted in. The first type in the list will be ordered before the second, 
				the second type before the third, etc.
				
Output:
-'keyList': a list of keys sorted in the order that the lines will be sorted.
-'holder': a dictionary of lines paired with keys in 'keyList'.
"""		
def sortGff3(gff3_File, types = ['gene','mRNA','exon']):
	f1 = open(gff3_File, "r")
	holder = dict()
	counter = dict()
	contigCount = False
	lineCount = -1
	formatCount = 0
    
	for line in f1:
        
		lineCount += 1
        
		try:
			theLine = line.strip().split("\t")
			if len(theLine) > 9:
				raise FormatError("0200")
			elif len(theLine) < 9:
				raise FormatError("0100")
		except FormatError as er:
			Errors.append("[" + str(lineCount) + "] " + er.returnError())
			continue
		except:
			Errors.append("Validation Error: unkown error. Check line: " + line)
        
		if theLine[2] == "contig":
			continue
        
		if theLine[2] not in types:
			Errors.append("The third component of each line must be one of the types. The types are: " + str(types) + " .")
			return "kill"
        
		holder[theLine[2]+"_"+theLine[3]] = line
        
		if theLine[3] in counter:
			counter[theLine[3]] = counter[theLine[3]] + types.index(theLine[2])
		else:
			counter[theLine[3]] = types.index(theLine[2])
			
	f1.close()
    
	counterKeys = counter.keys()
    
	expectedNum = ((len(types))*(len(types)-1))/2
	for key in counterKeys:
		thisNum = counter[key]
		try:
			if thisNum == expectedNum:
				continue
			else:
				raise FormatError("0500")
		except FormatError as er:
			Errors.append("Coordinate " + key + " " + er.returnError())
    
	keyList = holder.keys()  
	i = 0    
	length = len(keyList)
	while i < length:
		Key = keyList[i].split("_")
      
		try:
			if Key[0] == "contig":
				if not contigCount:
					temp = keyList[0]
					keyList[0] = keyList[i]
					keyList[i] = temp
					contigCount = True
					continue
				else:
					raise FormatError("0400") ## should only be one contig per file
		except FormatError as er:
			Errors.append("[" + str(i) + "] " + er.returnError())
			i+=1
			continue ## may need to look at this closer. Handling when there are multiple contig lines.
        
		if contigCount:
			if i == 1:
				i+=1
		else:
			if i == 0:
				i+=1
            
		preKey = keyList[i-1].split("_")
        
		if preKey[1] == Key[1]: ## same start coordinates
			prePos = int(types.index(preKey[0]))
			Pos = int(types.index(Key[0]))
			if Pos < prePos : ## sort based on types list. Start of list > end.
				temp = keyList[i-1]
				keyList[i-1] = keyList[i]
				keyList[i] = temp
				i-=1
				continue
			else:
				i+=1
				continue
        
		elif int(preKey[1]) > int(Key[1]): ## preveous key's coords is > current key's coord
			temp = keyList[i-1]
			keyList[i-1] = keyList[i]
			keyList[i] = temp
			i-=1
            
		else:
			i+=1
            
	return [keyList, holder]
    
"""
Checks each component of each line of the gff file for proper format. Prints to the global
list 'Errors' when a component is incorrectly formatted.

Parameters:
-'keyList': list of sorted keys.
-'holder': dictionary of lines paired with keys in 'keyList'.
-'Seq': string of the nucleotide sequence extracted from uploaded fasta file.
-'types': a list indicating the types of each gff line and the order the types should
				be sorted in. The first type in the list will be ordered before the second, 
				the second type before the third, etc.

"""    
def fileCheck(keyList, holder, Seq, types = ['gene','mRNA','exon']):

	priorID = "N/A"
	count = 0
	lineCount = 0
	namesList=[]
	
	for key in keyList:
	    
		lineCount += 1
        
		try:
			theLine = holder[key].strip().split("\t") ## try here to catch if students not tab deliminating or adding extra lines
		except:
			Errors.append("Format Error: each line needs to be tab deliminated.")
			continue
        
		try:
			if len(theLine) > 9:
				raise FormatError("0200")
			elif len(theLine) < 9:
				raise FormatError("0100")
		except FormatError as er:
			Errors.append("[" + str(lineCount) + "] " + er.returnError())
			continue
            
        ### 1ST and 2ND ITEM ### 
		try:
			if charCheck(theLine[0]):
				raise LineError("0002")
			elif charCheck(theLine[1]):
				raise LineError("0003")
		except LineError as er:
			Errors.append("[" + str(lineCount) +"] " + er.returnError())
        
        ### 3RD ITEM ###
		try:
			if theLine[2] not in types:
				raise LineError("0004")  ## types can be changed if more types of line needed
		except LineError as er:
			Errors.append("[" + str(lineCount) +"] " + er.returnError())
        
        ### 4TH ITEM ###
		try:
			if not int(theLine[3]) > 0:
				raise LineError("0005") ## needs to be greater than 0
		except LineError as er:
			Errors.append("[" + str(lineCount) +"] " + er.returnError())
		except ValueError as er:
			Errors.append("[" + str(lineCount) +"] " + str(er)) 
        
        
        ### 5TH ITEM ###
		try:
			if not int(theLine[4]) > int(theLine[3]):
				raise LineError("0006") ## needs to be greater than 1st coordinate
		except LineError as er:
			Errors.append("[" + str(lineCount) +"] " + er.returnError())
		except ValueError as er:
			Errors.append("[" + str(lineCount) +"] " + str(er)) 
        
        ### 6TH ITEM ###
		try:
			if theLine[5] != ".": ## can also be a number
				try:
					int(theLine[5])
				except ValueError:
					raise LineError("0007")
		except LineError as er:
			Errors.append("[" + str(lineCount) +"] " + er.returnError())
            
        ### 7TH ITEM ###
		try:
			if theLine[6] != "+":
				if theLine[6] != "-":
					raise LineError("0008")  ## + or - strand
		except LineError as er:
			Errors.append("[" + str(lineCount) +"] " + er.returnError())
                
         ### 8TH ITEM ### 
		try:
			if theLine[7] != ".":
				raise LineError("0009")
		except LineError as er:
			Errors.append("[" + str(lineCount) +"] " + er.returnError())
                
        
        ### 9TH ITEM ###        
		if theLine[2] == types[0]: ######################## FOR types[0] (Where the gene is checked) #############
			geneCheck(int(theLine[3]),int(theLine[4]),Seq,lineCount)
			last = theLine[8]
            
			try:
				if charCheck(last) :
					raise LineError("0011")
			except LineError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError())
            
			count=1
            
			last = theLine[8].split(";")
            
			_Name = False
			_ID = False
			try:
				for x in range(0, len(last)):
					if last[x].find("=") < 0:
						continue
					elif last[x].find("ID=") == 0:
						_ID = True
						priorID = last[x][3:]
						if namesList.count(priorID) > 0:
							raise LineError("0012") ## each id can only be used once 
						else:
							namesList.append(priorID)
					elif last[x].find("Name=") == 0:
						_Name = True
					else:
						raise ValidationError("1000")
						continue
			except LineError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError())
			except ValidationError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError() + " = spurious info. Recheck requirements of 9th component")
            
			try:
				if not _Name or not _ID:
					raise LineError("0013")  ## has to have a Name and an ID
			except LineError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError())
                
                
                                         ###################### FOR ALL OTHER TYPES ##################
		else:                      
			typesPos = int(types.index(theLine[2]))
            
			try:
				if not count == typesPos:
					raise LineError("0021") ## Either the file is not sorted properly or not all types are present for each gene
					continue
			except LineError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError())
                
			count += 1
			last = theLine[8]
            
			try:
				if charCheck(last):
					raise LineError("0022")
			except LineError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError())
            
			tempID = "N/A"
			last = theLine[8].split(";")
            
			_Parent = False
			_ID = False
			try:
				for x in range(0, len(last)):
					if last[x].find("=") < 0:
						continue
					elif last[x].find("ID=") == 0:
						_ID = True
						tempID = last[x][3:]
						if namesList.count(tempID) > 0:
							raise LineError("0023") ## each id can only be used once 
						else:
							namesList.append(tempID)
					elif last[x].find("Parent=") == 0:
						_Parent = True
					else:
						raise ValidationError("1000") ## Catch everything else.
			except LineError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError())
			except ValidationError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError() + " = spurious info. Recheck requirements of 9th component")
				continue
            
			try:
				if count == len(types):
					if not _Parent:
						raise LineError("0024") ## for last type must at least have a Parent
				else:
					if not _Parent or not _ID:
						raise LineError("0025") ## need an ID and Parent
			except LineError as er:
				Errors.append("[" + str(lineCount) +"] " + er.returnError())
                
	return "clean"
    
"""
Checks a string for characters not a-zA-Z0-9.=; and returns true if such a character is
found.

Parameters:
-'str': string being checked.
"""
def charCheck(str, search=re.compile(r'[^a-zA-Z0-9.=;]').search):
	return bool(search(str))
    
"""
Checks if coordinates of sequence given is a proper gene. Writes errors to global list
'Errors'.

Proper gene:
- can be translated to a protein.
- has a stop codon.
- no internal stop codon.
- begins with a start codon.

Parameters:
-'coord1': start coordinate of gene.
-'coord2': second coordinate of gene.
-'Seq': string nucleotide sequence.
-'count': line in the file the gene is from.
"""
def geneCheck(coord1, coord2, Seq, count):
	gene = Seq[coord1-1:coord2]
	try:
		if len(gene)%3 != 0:
			raise BiologyError("0050") ## divisible by three
	except BiologyError as er:
		Errors.append("[" + str(count) +"] " + er.returnError())
    
	try:
		protein = translate(gene) ## errors in translating the gene
	except:
		Errors.append("[" + str(count) + "]  Biology Error: unable to translate sequence. Ensure sequence provided is a nucleotide sequence.")
		return

	Start = True
	startCodon = Seq[coord1-1:coord1+2]
	if startCodon == "ATG" or startCodon == "TTG" or startCodon == "GTG":
		Start = False
    
	try:
		if protein.count("*") == 0:
			raise BiologyError("0030") ## no stop codon
		if protein.count("*") > 1:
			raise BiologyError("0040") ## internal stop codons
		if Start:
			raise BiologyError("0020") ## has to start with M
	except BiologyError as er:
		Errors.append("[" + str(count) +"] " + er.returnError())
    
	return
    
"""
Reads in a nucleotide sequence from a text file in .fasta format.

Parameters:
-'fasta_File': uploaded text file in .fasta format

Output:
-'seq': string of a nucleotide sequence.
"""
def fastaRead(fasta_File):

	seq = []
	try:
		f1 = open(fasta_File, "r")
	except:
		return "Error: unable to read fasta file."
    
	#for line in fasta_File:
	for line in f1:
	
		line = line.rstrip()
		if line.startswith(">"):
			continue
		else:
			seq.append(line)

	f1.close()
	return ''.join(seq)

"""
Translates a nucleotide sequence into a protein sequence. Uses a dictionary of nucleotides
paired with the protein they translate for.

Parameters:
-'nucSeq': string of nucleotides.

Output:
-'protSeq': string of proteins.
"""
def translate(nucSeq):
    
	codonLib = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V',
	'GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A',
	'GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K',
	'GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R',
	'AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    
	cntLoc = 0
	cntCodon = 0
	codon = ''
	protSeq = ''
	while cntLoc < len(nucSeq):
		codon = codon + nucSeq[cntLoc]
		cntCodon += 1
		if cntCodon == 3:
			temp = codonLib[codon]
			protSeq = protSeq + temp
			codon = ''
			cntCodon = 0
		cntLoc += 1
        
	return protSeq

"""
Writes a sorted gff file and an errors text file.

Parameters:
-'keyS': list of sorted keys.
-'holderS': dictionary of lines from the gff file paired with keys in 'keyS'.
-'nameS': name of new sorted gff file. Also location where the file is written.
-'errors': list of all errors found in uploaded files to be written to text file.
-'nameE': name of new errors text file. Also location where the file is written.
-'incLine': boolean indicating whether lines from the sorted file should be included in 
				errors file.
"""
def outFile(keyS, holderS, nameS, errors, nameE, incLine):

	outE = nameE
	outS = nameS
	
	if incLine:
		for item in errors:
			if "Coordinate" in item:
				tempE = item.split(" ")[1]
				for key in keyS:
					if tempE == key.split("_")[1]:
						outE.write(holderS[key])
						break
			else:
				tempE = int(item[1])
				outE.write(holderS[keyS[tempE]])
			
			outE.write(item)
			outE.write("\n")
			outE.write("\n")
	else:
		for item in errors:
			outE.write(item)
			outE.write("\n")
	
	for key in keyS:
		outS.write(holderS[key])
		outS.write("\n")		
		


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
    

if __name__ == "__main__":
	main()
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    