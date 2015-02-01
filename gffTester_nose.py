import re,sys

def validHeader (line):
    """
    Return True if line is a valid gff3 header.
    
    Parameters:
    - 'line' - line to check
    """
    if line == "##gff-version 3\n":
        return True
    else:
        return False
    
def validTabStructure(line):
    """
    Checks a line for valid <tab> structure, a valid GFF3 file should have 9 entries and 8 <Tab> charachters
    
    Parameters:
    - 'line': line to check
    """
    theLine = line.strip().split("\t")  #strip off extra charachters and split the line on tab
    if len(theLine) != 9:               # 9 columns in a gff file, so return False if not exactly 9
        return False
    else:
        return True

def validSeqname(text):
    """
    Checks a string to make sure it is a valid entry for column 1 of gff file, for this column it must match
    the sequence name on the Gbrowse database and have only valid characters
    
    Parameters:
    - 'text': string value of entry in column 1 to check
    """
    return not charCheck(text)  # need name check, for now just check for valid charachters

def validSource(text):
    """
    Checks a string to make sure it is a valid entry for column 2 of gff file, for this column should be source
    which only need validity of the characters
    
    Parameters:
    - 'text': string value of the entry in column 2 to check
    """
    return not charCheck(text)

def validType(text, validTypes=None):
    """
    Checks a string to make sure it is a valid entry for column 3 of gff file, for this column should be feature type.
    For phage this should be one of types ['gene','mRNA','exon']
    
    Parameters:
    - 'text': string value of the entry in column 2 to check
    - 'validTypes': List of valid types for checking if not the default 3
    """
    if validTypes is None:
        validTypes = {'gene','mRNA','exon'}
    
    return text in validTypes


def validScore(score):
    """
    Checks a score (entries in column 6), should be a number or "."
    
    Parameters:
    - 'score': score from column 6
    """
    if score == ".":
        return True
    
    try:
        x = float(score)
    except:
        isScore = False
        pass
    else:
        isScore = True        

    return isScore

def validStrand(strand):
    """
    Checks a strand (entries in column 7), should be one of ("=", "-", ".", "?")
    
    Parameters:
    - 'strand': strand entry from column 7
    """
    validValues = {"+", "-", ".", "?"}
    
    return strand in validValues

def validPhase(phase):
    """
    Checks a phase (entries in column 8), should be one of (".", 0, 1, 2, ")
    
    Parameters:
    - 'phase': strand entry from column 8
    """
    validValues = {".", 0, "0" , 1, "1", 2, "2"}
    
    return phase in validValues

def validAttributes(attributes):
    """
    Checks a column 9 entry, should be series of key=value entries separated by ; 
    It is OK to have spaces in the values entry but not the Key value
    
    
    Parameters:
    - 'attributes': entire entry from column 9
    """
    
    #ok to have a null string
    if attributes == '':
        return True
    
    #checking for invalid characters goes here if any found
    
    if bool(charCheck(attributes, search=re.compile(r'[^a-zA-Z0-9.=;_ ]').search)):
        return False
    
    
    #ok go ahead and split into the underlyine key=value items 

    attrPairs = attributes.split(";")   #list of each attribute pair
    
    validity = True
    
    for attrPair in attrPairs:
        if attrPair == '':                  #this will happen if file had two ;; in a row
            continue
            
        if attrPair.count('=') != 1:        #for each Key=value there must be only one "="
            return False
        
        attrKey, attrValue = attrPair.split('=')
        
        if len(attrKey) < 1:                #must have an entry for a key
            validity = False
        
        if charCheck(attrKey):              #check key for valid charachters
            validity = False
        
        # check for invalid charachters in value, but spaces are OK, so remove before checking
        if charCheck(attrValue.replace(' ','')):
            validity = False
            
    return validity

def charCheck(str, search=re.compile(r'[^a-zA-Z0-9.=;_]').search):
    """
    Checks a string for characters NOT a-zA-Z0-9.=; and returns True if invalid character is
    found.

    Parameters:
    -'str': string being checked.
    -'search': default search all characters in valid set
    """
    return bool(search(str))

def validCoordinate(coord):
    """
    Checks a coordinate (entries in column 4 or 5), should be a positive integer
    
    Parameters:
    - 'coord': coordinate from column 4 or 5
    """
    
    return not bool(charCheck(coord, search=re.compile(r'[^0-9]').search))

def validCoordinates(leftCoord, rightCoord):
    """
    Checks a coordinates (entries in column 4 and 5), should be a positive integers
    and the left Coordinate should be smaller than the right Coordinate
    
    Parameters:
    - 'leftCoord':  coordinate from column 4
    - 'rightCoord': coordinate from column 5

    """
    return (validCoordinate(leftCoord) and validCoordinate(rightCoord) and int(leftCoord) <= int(rightCoord))


def printFailureMessage(failType):
    print "##### Fatal Error #####"
    print failType + "failed. No other tests run."
    print """
    
Computer messages follow below OK to ignore, you should fix file and try again.
================================ End Results ================================    
    """