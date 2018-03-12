import numpy as np


class Molecule:
    """ Class doc """
    
    def __init__ (self, index = None, name = 'ligand',  atoms = [], bonds = [], sectors = [], connections = []):
        """ Class initialiser """
        self.index     = index
        self.atoms     = atoms
        self.name      = name
        
        
        
        self.bonds     = bonds
        self.angles    = []
        self.dihedrals = []
        
        self.sectors   = []
        self.connections = connections
        
	
	# atributos utilizados no codigo antigo
        #--------------------------------------------------------------
	
	#                                         (int)            (a list of ints)                                 
	self.graph                   = {}  # {atom_x.index:  [atom_y.index, atom_z.index] , bla bla bla}
        

	#                                        (int , int) tupla        (a bond obj)                                 
	self.graph_bonds             = {}  # {(atom_x.index, atom_y.index) :  bond }
	
	#                   
	self.rotatable_bonds_edited  = []
	self.block_connections = {}
	#--------------------------------------------------------------
	
    def building_graph (self):
        """ Function doc """
        
        # construindo o grafo
	#--------------------------------------------------------------
	for atom in self.atoms:
            #print atom.bonds
            self.graph[atom.index] = []

	    #atribiur o proprio objeto "molecula" ao atributo "molecula" do objeto "bond"
	    atom.molecule = self
	    
	    
	for bond in self.bonds:
            self.atoms[bond.atom1.index].bonds.append(bond)
            self.atoms[bond.atom2.index].bonds.append(bond)
            self.graph[bond.atom1.index].append(bond.atom2.index)
            self.graph[bond.atom2.index].append(bond.atom1.index)
	    #atribiur o proprio objeto "molecula" ao atributo "molecula" do objeto "bond"
	    bond.molecule = self
	    
	    # buildind self.graph_bonds  
	    self.graph_bonds[(bond.atom1.index, bond.atom2.index)] = bond
	    self.graph_bonds[(bond.atom2.index, bond.atom1.index)] = bond
	
	
	print '''\n\ncontrucao do grafo self.graph'''
        print self.graph
	
	print '''\n\ncontrucao do grafo self.graph_bonds'''	
	print  self.graph_bonds
	#--------------------------------------------------------------

        #------------------------------------------------------------------------------------------------------------#
        # finding atoms with only one connection                                                                               
	print '\n\n ------ Finding atoms with only one connection ------ '
	not_alone_atoms = []                                                   
        for atom in self.atoms:                                                 
            #degree = self.Graph.vertex_degree(key)                             
            print(atom.index, len(atom.bonds),atom.name )                                                 
            if len(atom.bonds) > 1:                                                     
                not_alone_atoms.append(atom)                                    
	print ' ------ -------------------------------------- ------ \n\n'
	#------------------------------------------------------------------------------------------------------------#
        
        
        #------------------------------------------------------------------------------------------------------------#
        rotatable_bonds = []
        _buffer = []                                                           
        for atom1 in not_alone_atoms:                                            
            _buffer.append(atom1)                                                
            
            for atom2 in not_alone_atoms:                                       
                if atom2 in _buffer:                                            
                    pass                                                       
                else:                                                          
                    path = self.find_all_paths(atom1.index, atom2.index)                
                    if len(path) == 1:                                         
                        if len(path [0]) == 2:                                 
                            rotatable_bonds.append([atom1.index, atom2.index])
        
        
	print '\n\n ------  rotatable_bonds    ------ '
	#print rotatable_bonds
	
	for rotatable_bond in rotatable_bonds:
	    atom1 = rotatable_bond[0]
	    atom2 = rotatable_bond[1]

	    print atom1, self.graph_bonds[(atom1, atom2)].atom1.name," ---- ", atom2, self.graph_bonds[(atom1, atom2)].atom2.name, ' -> bond index ',  self.graph_bonds[(atom1, atom2)].index, self.graph_bonds[(atom1, atom2)].bond_type 
	
	print ' ------ --------------- ------ \n\n'

	#------------------------------------------------------------------------------------------------------------#
    
    
    
    
        # --------------------------------------------------------
        '''
        finding the methyl groups
        '''
        _buffer = []
        
        methyl_groups = []
        for bond in rotatable_bonds:
            for vertex_v in bond:
                if vertex_v in _buffer:
                    pass
                else:
                    _buffer.append(vertex_v)
                    n = 0 
                    
                    for vertex_w in self.graph[vertex_v]:
                        
                        if self.atoms[vertex_w].name == 'H':
                            #print vertex_v, vertex_w, 'H'
                            n += 1
                    
                    if n == 3:
                        print vertex_v, self.atoms[vertex_v].name, 'methyl group'
                        methyl_groups.append(vertex_v)
        # --------------------------------------------------------
        
        
        
        # --------------------------------------------------------
        print '''excluding methyl groups from rotatable bonds'''
	n = 0
        for bond in rotatable_bonds:
            if bond[0] in methyl_groups or bond[1] in methyl_groups:
                n+=1
		pass
            else:
                self.rotatable_bonds_edited.append(bond)
        # --------------------------------------------------------
	print 'total  numer of excluded methyl groups:', n
    
    

        print ('''
#-------------------------------------------
Looking for rotatble block / sectors        
#-------------------------------------------
	'''
	
	)
	self.find_rotatable_blocks2(rot_bonds = self.rotatable_bonds_edited)
        for block in self.blocks:
            self.block_connections[block] = []
        
        
        for bond in self.rotatable_bonds_edited:
            
            blocka = None
            blockb = None 
            
            for block in self.blocks:
                if bond[0] in self.blocks[block]:
                    #print bond[0], block,  self.Graph.blocks[block]
                    blocka = block
                if bond[1] in self.blocks[block]:
                    #print bond[1], block,  self.Graph.blocks[block]
                    blockb = block

            self.block_connections[blocka].append([blockb,[bond[0], bond[1]]])
            self.block_connections[blockb].append([blocka,[bond[1], bond[0]]])
       
	
	#-------------------------------------------------------------------------------------------
	print '\nNumber of sectors(blocks): ', len(self.blocks)
	#print self.blocks
	
	for sector_index in self.blocks:
	    #print sector_index, self.blocks[sector_index]
	    
	    
	    sector = Sector(index = sector_index -1, 
	                    label = 'Sector', 
			    atoms = [], 
			    bonds = [], 
			    connections = [], 
			    molecule = self
			    ) 
	    for atom_index in self.blocks[sector_index]:
		sector.atoms.append(self.atoms[atom_index])
	    
	    for atom in sector.atoms:
		print sector.index, '---> ',atom.index, atom.name, atom.symbol, atom.atom_type, atom.pos, len(atom.bonds) #, atom.charge, atom.sector.index, len(atom.bonds), atom.molecule.name 
		atom.sector = sector
		
	    self.sectors.append(sector)
	    
	    
	print '\nConnections of sectors(blocks): ', len(self.blocks)
	
	
	
	for connection_index1 in self.block_connections:
	    
	    #print connection_index1, self.block_connections[connection_index1]
	    
	    for connection_data in  self.block_connections[connection_index1]:
		print connection_index1, connection_data[0], connection_data[1]
		
		#'''
		connection =  Connection(sector1  = self.sectors[connection_index1 -1]  , 
					 sector2  = self.sectors[connection_data[0]-1] , 
					 bond     = self.graph_bonds[(connection_data[1][0],connection_data[1][1])], 
					 molecule = self
					 )	#'''
	
		
		self.sectors[connection_index1 -1].connections.append(connection)
		#self.sectors[connection_data[0]-1].connections.append(connection)
		self.connections.append(connection)
	
	for sector in self.sectors:
	    
	    for connection in  sector.connections:
		print 'sector', sector.index ,'sector a:', connection.sector1.index, 'sector b:',connection.sector2.index, 'atom a:',connection.bond.atom1.index,connection.bond.atom1.name,'atom b:' ,connection.bond.atom2.index ,connection.bond.atom2.name,'bond index:',connection.bond.index
	
	#-------------------------------------------------------------------------------------------
	
	
	print 'apague isso depois'
	for atom in self.atoms:
	    print atom.sector.index, len(atom.sector.atoms)
	
	#for  sector in  easydock.ligands[0].sectors:
	#    for atom in sector.atoms:
	#	print sector.index, atom.index, atom.name
	#	atom.sector = sector 	

    def find_rotatable_blocks2 (self, rot_bonds = []):
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
	
	blocks = {}
	graph  = self.graph
	
	edited_graph = graph
	_buffer      = [] 
	
	for bond in rot_bonds:
	    
	    edited_graph[bond[0]].remove(bond[1])
	    edited_graph[bond[1]].remove(bond[0])
	    if bond[0] in _buffer :
		pass
	    else:
		_buffer.append(bond[0])
		
	    if bond[1] in _buffer :
		pass
	    else:
		_buffer.append(bond[1])
	
	#print graph
	#print edited_graph
	#print _buffer
	
	
	#print vertex_v
	#for key in self.visited:
	#    self.visited[key] = False
	
	
	#self.blocks  = {}
	self.block_index = 1
	self.blocks  = {}
	
	
	for vertex_v in _buffer:
	    
	    #'''
	    vertex_v_already_in_blocks =  False
	    #'''
	    
	    #'''
	    #---------------------------------------------------
	    # verifica se o vertex jah faz parte de algum bloco
	    #---------------------------------------------------
	    for block_index  in self.blocks:
		 
		 #print 'vertex:',vertex_v, 'block:', block_index, 'items', self.blocks[block_index]
		 
		 if vertex_v in self.blocks[block_index]:
		     #print vertex_v, self.blocks[self.block_index]
		     vertex_v_already_in_blocks = True
	    #---------------------------------------------------
	    #'''
	    
	    '''
	    #print vertex_v_already_in_blocks
	    print self.blocks
	    for block_index  in self.blocks:
		print 'block:', block_index, 'items', self.blocks[block_index]
	    #'''
	    
	    if vertex_v_already_in_blocks:
		pass
	    
	    else:		
		self.blocks[self.block_index] = []
		#print 'vertex', vertex_v
		self.busca_profunda ( edited_graph, vertex_v)
		#print '\n'
	    
		self.block_index += 1
		self.blocks[self.block_index] = []
		 
    def busca_profunda (self, grafo, vertex_v):
	""" Function doc """
	
	#block  = []
	#block.append(vertex_v)
	#print vertex_v
	
	if vertex_v in self.blocks[self.block_index]:
	    pass
	else:
	    self.blocks[self.block_index].append(vertex_v)
	
	#for key in self.visited:
        #    self.visited[key] = False
	
	for atom in self.atoms:
	    atom.visited = False
	
	
	
	for vertex_w in grafo[vertex_v]:
	    if self.atoms[vertex_w].visited == False:
		self.busca_profunda2(grafo, vertex_w)

    def busca_profunda2 (self, grafo, vertex_w):
	""" Function doc """

        self.atoms[vertex_w].visited = True
	#print vertex_w
	
	
	if vertex_w in self.blocks[self.block_index]:
	    pass
	else:
	    self.blocks[self.block_index].append(vertex_w)
	
	for vertex_k in grafo[vertex_w]:
		
	    if self.atoms[vertex_k].visited == False:
                self.busca_profunda2(grafo, vertex_k)
    
    def printAtoms (self):
        """ Function doc """
        for atom in self.atoms:
            print atom.index, atom.name, atom.symbol, atom.atom_type, atom.pos, len(atom.bonds) #, atom.charge, atom.sector.index, len(atom.bonds), atom.molecule.name 
        pass

    def printBonds (self):
        """ Function doc """
        for bond in self.bonds:
            print bond.atom1.index, bond.atom2.index, bond.bond_type

    def find_path(self, start_vertex, end_vertex, path=None):
        """ find a path from start_vertex to end_vertex 
            in graph """
        if path == None:
            path = []
        
        graph = self.graph
        
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
        graph = self.graph
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








class Connection:
    """ Class doc """
    
    def __init__ (self, sector1 = None, sector2 = None, bond = None, molecule = None):
        """ Class initialiser """
        self.sector1  = sector1
        self.sector2  = sector2
        self.bond     =  bond
        self.molecule = molecule
        self._type    = 'simple'
        
        
class Sector:
    """ Class doc """
    
    def __init__ (self, index = None, label = 'Sector', atoms = [], bonds = [], connections = [], molecule = None):
        """ Class initialiser """
        self.index       = index
        self.label       = 'Sector'
        self.actived     = False
        self.atoms       = atoms
        
        self.bonds       = bonds
        self.connections = connections
        self.molecule    = molecule
        pass


class Bond:
    """ Class doc """
    
    def __init__ (self, index = None, atom1 = None, atom2 = None, b_type = 1, sector = None, molecule = None):
        """ Class initialiser """
        #self.atoms     = []
        self.index     = index
        self.atom1     = atom1
        self.atom2     = atom2
        self.bond_type = b_type
        self.sector    = sector
        self.molecule  = molecule
        

class Atom:
    """ Class doc """
    
    def __init__ (self, index        = None, 
                        name         = 'Xx',
                        symbol       = 'X' , 
                        atom_type    = None, 
                        pos          = None, 
                        charge       = 0.00,

                        sector       = None,
                        resi         = None, 
                        resn         = None, 
                        chain        = None, 
                        
                        bonds        = [] ,
                        molecule     = None,
                        ):
 
        if pos is None:
            pos = np.array([0.0, 0.0, 0.0])

        self.index      = index     
        self.name       = name      
        self.symbol     = symbol    
        self.atom_type  = atom_type 
        self.pos        = pos       
        self.charge     = charge    
                    
        self.sector     = sector    
        self.resi       = resi      
        self.resn       = resn      
        self.chain      = chain     
                        
        self.bonds      = bonds 
        self.molecule   = molecule  
        
        self.visited    = False
        
        
        
        
        
        
