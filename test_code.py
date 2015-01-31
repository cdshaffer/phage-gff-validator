
from gffTester_nose import *

'''
need to check the following functions:

def validHeader (line): TEST valid, extra space, extra tab, garbled
def validTabStructure(line): test: 7 tabs, 8 tabs, 9 tabs
def validSeqname(text): SKIP, just a "not charCheck" for now
def validSource(text): SKIP, just a "not charCheck" for now
def validType(text, validTypes=None): TEST: pass with any 3 defaults 'gene','mRNA','exon'; fails on none item, passes with list
def validCoordinate(coord): TEST: fails with text, float, negative int; passes with positive int
def validScore(score):
def validStrand(strand):
def validPhase(phase):
def validAttributes(attributes):
def charCheck(str, search=re.compile(r'[^a-zA-Z0-9.=;_]').search):
def validGene(seq, line):       
def fastaRead(fasta_File):
def countStopCodons(seq, strand):
'''

def test_validHeader_1():
    'valid header tests as valid'
    result = validHeader('##gff-version 3\n')
    assert result == True

def test_validHeader_2():
    'header with extra space fails'
    result = validHeader('##gff-version 3 \n')
    assert result == False
    
def test_validHeader_3():
    'header with extra tab fails'
    result = validHeader('##gff-version 3\t\n')
    assert result == False

    
def test_validHeader_4():
    'garbled header fails'
    result = validHeader('##gff-Version 3\n')
    assert result == False
    
def test_validTabStructure_1():
    'fails 7 tab structure'
    result = validTabStructure("1\t2\t3\t4\t5\t6\t7\t8\n")
    assert result == False

def test_validTabStructure_2():
    'passes 8 tab structure'
    result = validTabStructure("1\t2\t3\t4\t5\t6\t7\t8\t9\n")
    assert result == True
    
def test_validTabStructure_3():
    'fails 9 tab structure'
    result = validTabStructure("1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n")
    assert result == False
    
def test_validType_1():
    'passes with "gene"'
    result = validType('gene')
    assert result == True

def test_validType_2():
    'passes with "mRNA"'
    result = validType('mRNA')
    assert result == True

def test_validType_3():
    'passes with "exon"'
    result = validType('exon')
    assert result == True
    
def test_validType_4():
    'fails with mrna'
    result = validType('mrna')
    assert result == False
    
def test_validType_5():
    'passes with specified list'
    result = validType('mrna', ['mRNA', 'mrna'])
    assert result == True
    

def test_validCoordinate_1():
    'fails with text as coordinate'
    result = validCoordinate('a')
    assert result == False

def test_validCoordinate_2():
    'fails with float'
    result = validCoordinate('2.0')
    assert result == False
    
def test_validCoordinate_3():
    'fails with negative int'
    result = validCoordinate('-3')
    assert result == False
    
def test_validCoordinate_4():
    'passes with positive integer'
    result = validCoordinate('100')
    assert result == True
    
def test_validCoordinates_1():
    'fails with left larger than right'
    result = validCoordinates('5','4')
    assert result == False

def test_validCoordinates_2():
    'passes with left equal to right'
    result = validCoordinates('5','5')
    assert result == True
    
def test_validCoordinates_3():
    'passes with left smaller than right'
    result = validCoordinates('4','5')
    assert result == True