from collections import defaultdict
import random

class DeBruijnGraph():
    """Main class for De Bruijn graphs
    
    Private Attributes:
        graph (defaultdict of lists): Edges for De Bruijn graph
        first_node (str): starting position for traversing the graph
    """

    def __init__(self, input_string, k):
        self.graph = defaultdict(list)
        self.first_node = ''
        self.build_debruijn_graph(input_string, k)
        
    def add_edge(self, left, right):
        ''' This function adds a new edge to the graph
        
        Args:
            left (str): The k-1 mer for the left edge
            right (str): The k-1 mer for the right edge

        Updates graph attribute to add right to the list named left in defaultdict   
        '''
        self.graph[left].append(right)
        
    def remove_edge(self, left, right):
        ''' This function removed an edge from the graph
        
        Args:
            left (str): The k-1 mer for the left edge
            right (str): The k-1 mer for the right edge

        Updates graph attribute to remove right from the list named left in defaultdict
        '''
        matching_edges = []
        for i, key in enumerate(self.graph[left]):
            if key == right:
                self.graph[left].pop(i)
                break

        
    def build_debruijn_graph(self, input_string, k):
        ''' This function builds a De Buijn graph from a string
        
        Args:
            input_string (str): string to use for building the graph
            k (int): k-mer length for graph construction

        Updates graph attribute to add all valid edges from the string
        
        Example:
        >>> dbg = DeBruijnGraph("this this this is a test", 4)
        >>> print(dbg.graph) #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        defaultdict(<class 'list'>, {'thi': ['his', 'his', 'his'], 'his': ['is ', 'is ', 'is '], ...)
        '''
        for i in range(len(input_string) - k + 1):
            kmer = input_string[i:i+k]
            left_mer = kmer[0:k-1]
            right_mer = kmer[1:k]
            self.add_edge(left_mer, right_mer)
            
            if i == 0:
                self.first_node = left_mer
        
        #self.add_edge(right_mer, self.first_node)
            
    def print_eulerian_walk(self, seed=None):
        ''' This function starts the recursive walk function
        at the first node in the graph
        
        Args: None
        
        Returns:
            tour (list): list of k-1 mers traversed by the algorithm
        
        Example:
        >>> dbg = DeBruijnGraph("this this this is a test", 4)
        >>> dbg.print_eulerian_walk(seed=1) #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        ['thi', 'his', 'is ', 's i', ' is', 'is ', ...]
        '''
        tour = []
        tour = self.eulerian_walk(self.first_node, seed=seed)
        tour.append(self.first_node)
        return tour[::-1]
        
    def eulerian_walk(self, node, seed=None):
        ''' This is a recursive function that follows all edges from a node
        to traverse the graph
        
        Args: 
            node (str): current node to traverse from
            seed (int): seed for random selection of edge to follow
        
        Returns:
            tour (list): list of k-1 mers traversed so far by the algorithm
            Note: this will be reverse order because of recursion
            
        Example:
        >>> dbg = DeBruijnGraph("this this this is a test", 4)
        >>> dbg.eulerian_walk('thi', seed=1) #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        ['is ', 'his', 'thi', ' th', ...]
        '''
        tour = []
        random.seed(seed)
        random.shuffle(self.graph[node])
        for next_node in self.graph[node]:
            self.remove_edge(node, next_node)
            tour = self.eulerian_walk(next_node, seed=seed)
            tour.append(next_node)
        return tour

    
#FASTA READER

def get_fasta(file):
    """Generator to lazily get all the fasta entries from a fasta file

    Args:
        fasta_file (str): /path/and/name/to.fasta

    Yields:
        header (str): header sequence of fasta entry (includes '>')
        seq (str): concatenated string sequence of the fasta entry
    """
    
    if '.gz' in file:
        file_object = gzip.open(file, 'rt')
    else:
        file_object = open(file, 'rt')
        
    name=''
    seq=''
    for line in file_object:
        
        # Capture the next header, report what we have, and update
        if line.startswith('>') and seq: #not first seq
            name = name[1:] #removes the carrot
            yield name, seq
            name=line.strip()
            seq=''
            
        # Get to the first header
        elif line.startswith('>'):  #first seq
            name=line.strip()
            
        # Just add sequence if it is the only thing there
        else:
            seq+=line.strip()
            
    # At the end, return the last entries
    if name and seq: #last seq
            name = name[1:]
            yield name, seq
