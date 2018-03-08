""" A Python Class
A simple Python graph class, demonstrating the essential 
facts and functionalities of graphs.
"""
#from MOL2Files import *

class Graph(object):

    def __init__(self, graph_dict=None):
        """ initializes a graph object 
            If no dictionary or None is given, 
            an empty dictionary will be used
        """
        if graph_dict == None:
            graph_dict = {}
        self.__graph_dict = graph_dict
    
        
        self.visited = {}
        for key in self.__graph_dict:
            self.visited[key]  = False
    
    
	self.pnum   = 0 
	self.prenum = {}
    
    
    
    
    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def add_edge(self, edge):
        """ assumes that edge is of type set, tuple or list; 
            between two vertices can be multiple edges! 
        """
        edge = set(edge)
        (vertex1, vertex2) = tuple(edge)
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]

    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if {neighbour, vertex} not in edges:
                    edges.append({vertex, neighbour})
        return edges

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res


    def find_path(self, start_vertex, end_vertex, path=None):
        """ find a path from start_vertex to end_vertex 
            in graph """
        if path == None:
            path = []
        
        graph = self.__graph_dict
        
        path = path + [start_vertex]
        
        if start_vertex == end_vertex:
            return path
        
        if start_vertex not in graph:
            return None
        
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_path = self.find_path(vertex, 
                                               end_vertex, 
                                               path)
                if extended_path: 
                    return extended_path
        return None



    def find_all_paths(self, start_vertex, end_vertex, path=[]):
        """ find all paths from start_vertex to 
            end_vertex in graph """
        graph = self.__graph_dict 
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return [path]
        if start_vertex not in graph:
            return []
        paths = []
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_paths = self.find_all_paths(vertex, 
                                                     end_vertex, 
                                                     path)
                for p in extended_paths: 
                    paths.append(p)
        return paths


    def vertex_degree(self, vertex):
        """ The degree of a vertex is the number of edges connecting
            it, i.e. the number of adjacent vertices. Loops are counted 
            double, i.e. every occurence of vertex in the list 
            of adjacent vertices. """ 
        adj_vertices =  self.__graph_dict[vertex]
        degree = len(adj_vertices) + adj_vertices.count(vertex)
        return degree


    def find_isolated_vertices(self):
        """ returns a list of isolated vertices. """
        graph = self.__graph_dict
        isolated = []
        for vertex in graph:
            print(isolated, vertex)
            if not graph[vertex]:
                isolated += [vertex]
        return isolated

    def delta(self):
        """ the minimum degree of the vertices """
        min = 100000000
        for vertex in self.__graph_dict:
            vertex_degree = self.vertex_degree(vertex)
            if vertex_degree < min:
                min = vertex_degree
        return min
        
    def Delta(self):
        """ the maximum degree of the vertices """
        max = 0
        for vertex in self.__graph_dict:
            vertex_degree = self.vertex_degree(vertex)
            if vertex_degree > max:
                max = vertex_degree
        return max



    @staticmethod
    def erdoes_gallai(dsequence):
        """ Checks if the condition of the Erdoes-Gallai inequality 
            is fullfilled 
            dsequence has to be a valid degree sequence
        """
        if sum(dsequence) % 2:
            # sum of sequence is odd
            return False
        for k in range(1,len(dsequence) + 1):
            left = sum(dsequence[:k])
            right =  k * (k-1) + sum([min(x,k) for x in dsequence[k:]])
            if left > right:
                return False
        return True
    
    def density(self):
        """ method to calculate the density of a graph """
        g = self.__graph_dict
        V = len(g.keys())
        E = len(self.edges())
        return 2.0 * E / (V *(V - 1))


    
    
    
    def find_rotatable_blocks (self, v1, v2):
	""" Function doc 
	
	  o---o          o---o
	 /     \        /     \
	o       o------o       o
	 \     /v1    v2\     /
	  o---o          o---o
			      \
			       o---o   
			      /     \  
			     o       o
			      \     /  
			       o---o   
			       
	|----------|------------------!
	
	   block 1       block 2
	   (list)        (list)
	
	Find all the "o" that is bonded to v2 but not with v1 
    
	"""
	
	graph  = self.__graph_dict
	self.block = []
	
	for key in self.visited:
	    self.visited[key] = False

	self.visited[v2]  = True
	#print self.visited
	#print v1, v2, self.visited[v2]

	for vertex in graph[v2]:
	    if vertex == v1:
		pass
	    
	    else:
		#print v1, v2, vertex, graph[v2]
		
		if self.visited[vertex] == False:
		    self.busca_block_prof(vertex)
		    
		    #block2 += block

	
	print v1, v2, self.block
	#print self.visited
	#return block2
    
    
    
    def busca_block_prof (self, vertex):
	""" Function doc """
        self.visited[vertex] = True
	self.block.append(vertex)
	
	#print vertex, self.__graph_dict[vertex]
	
	for w in self.__graph_dict[vertex]:
            
	    if self.visited[w] == False:
                self.busca_block_prof(w)
	
	#return block
    
    
    
    
    
    
    
    def  Busca (self, excluded = None):
	'''
	procedimento Busca(G: Grafo)
	pnum := 0
	Para Cada vertice v de G:
	Marque v como nao visitado
	Para Cada vertice v de G:
	Se v nao foi visitado:
	Busca-Prof(v)
	'''

	graph = self.__graph_dict




	for key in self.visited:
	    self.visited[key] = False

	for vertex in graph.keys():
	    if vertex == excluded:
		pass
	    else:
		if self.visited[vertex] == False:
		    self.Busca_Prof(vertex)
	
	print self.prenum
    
    
    def  Busca_Prof(self, v):
        """  
        procedimento Busca-Prof(v: vertice)
            Marque v como visitado
            pnum := pnum + 1
            prenum[v] := pnum
            Para Cada vertice w adjacente a v:
                Se w nao foi visitado:
                    Busca-Prof(w)       
        """
        self.visited[v] = True
        self.pnum += 1
        self.prenum[v] = self.pnum
        
	print v, self.__graph_dict[v]
        
	for w in self.__graph_dict[v]:
            if self.visited[v] == False:
                self.Busca_Prof(w)

    
    def find_path(self, start_vertex, end_vertex, path=None):
        """ find a path from start_vertex to end_vertex 
            in graph """
        if path == None:
            path = []
        
        graph = self.__graph_dict
        
        path = path + [start_vertex]
        
        if start_vertex == end_vertex:
            return path
        
        if start_vertex not in graph:
            return None
        
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_path = self.find_path(vertex, 
                                               end_vertex, 
                                               path)
                if extended_path: 
                    return extended_path
        return None



    


    






if __name__ == "__main__":
    '''
    g = { "a" : ["b", 'c'],         #1
	      "b" : ["d", 'e', "a"],    #2
          "c" : ["a", "d"],         #3
          "d" : ["b", "e", 'f'],    #4
          "e" : ["b","c"]     ,     #5
          "f" : ['d', 'g', 'h' ],   #6
	      "g" : ['f', 'h']    ,     #7
	      "h" : ['f','g']           #8
	}                           



    g = { 1 : [2, 3],         #1
	      2 : [4, 5, 1],    #2
          3 : [1, 4],         #3
          4 : [2, 5, 6],    #4
          5 : [2, 3]     ,     #5
          6 : [4, 7, 8],   #6
	      7 : [6, 8]    ,     #7
	      8 : [6, 7]           #8
	}  
    '''

    #atoms, _graph = load_mol2_files (infile = '/home/fernando/programs/EasyDock/mol2/com7.mol2', VMSession =  None, gridsize = 3)


    #graph = Graph(g)
    #graph =Graph(g)
    #print graph.visited
    #graph.Busca()
    
    
    #print graph.edges()
    #print("Vertices of graph:")
    #print(graph.vertices())
    #
    #print("Edges of graph:")
    #print(graph.edges())
    #
    #print("Add vertex:")
    #graph.add_vertex("z")
    #
    #print("Vertices of graph:")
    #print(graph.vertices())
    #
    #print("Add an edge:")
    #graph.add_edge({"a","z"})
    #
    #print("Vertices of graph:")
    #print(graph.vertices())
    #
    #print("Edges of graph:")
    #print(graph.edges())
    #
    #print('Adding an edge {"x","y"} with new vertices:')
    #graph.add_edge({"x","y"})
    #print("Vertices of graph:")
    #print(graph.vertices())
    #print("Edges of graph:")
    #print(graph.edges())
    #
    
    
    
    #print graph.find_path("d", "f", path=None)
    #
    #print('All paths from vertex "c" to vertex "c":')
    #path = graph.find_all_paths("d", "f")
    #print(path)
    #
    #
    #print(graph.vertex_degree("c"))
    
    
    '''
    not_alone_atoms = []
    
    for key in _graph:
        degree = graph.vertex_degree(key)
        print(key, degree)
        if degree > 1:
            not_alone_atoms.append(key)

    _buffer = []
    for key in not_alone_atoms:
        _buffer.append(key)
        for key2 in not_alone_atoms:
            if key2 in _buffer:
                pass
            else:
                path = graph.find_all_paths(key, key2)
                if len(path) == 1:
                    if len(path [0]) == 2:
                        print(key+1, key2+1, len(path),len(path [0]) )    
    '''
