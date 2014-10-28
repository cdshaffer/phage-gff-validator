
from gffTester_nose import *

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
    