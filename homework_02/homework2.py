#function to caculate BWT string
def BWT(string):
    ''' Function to calculate Burrows-Wheeler Transform for a given string.
    
    Args:
        string (str): Input string
    
    Returns:
        bwt_string (str): BWT string        
        
    Example:
        >>> BWT('googol')
        'lo$oogg'
        
    '''
    #Append '$' to the end of string
    string += '$'
    
    #generate table of circulated strings
    t = sorted(string[i:] + string[:i] for i in range(len(string)))
    #concatenate last symbols of circulated strings to generate BWT string
    bwt_string = ''.join([l[-1] for l in t])

    return bwt_string
