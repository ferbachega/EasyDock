import graphs_class as gclass


class Molecule:
    """ Class doc """
    
    def __init__ (self, 
    
                        name              = None ,
                        size              = None ,
                        _type             = None ,
                        charge            = None ,
                        info              = None ,
                        index             = None ,

                        atoms             = []   , 
                        bonds             = []   , 
                        
                        #bonds
                        index_bonds_pairs = None ,
                        bond_type         = None , 
                        
                        ):
        """ Class initialiser """
        self.name       = name
        self.size       = size  
        self._type      = _type 
        self.charge     = charge
        self.info       = info  
        self.index      = index

        
        
        
        
        
        self.atoms      = atoms
        for atom in self.atoms:
            atom.molecule =  self
        
        # -------------- Bonds --------------------
        self.bonds             = bonds
        self.bond_type         = bond_type
        self.index_bonds_pairs = index_bonds_pairs
        # -----------------------------------------

        #-----------------------------------------------------------------------------------------------------------#
        self.rotatable_bonds = []
        
        self.graph = {}                                                        
        for atom in self.atoms:                                                
            self.graph[atom.index] = []                                        
                                                                               
        for bond in self.index_bonds_pairs:                                    
            self.graph [ bond[0]].append(bond[1])                              
            self.graph [ bond[1]].append(bond[0])                              
                                                                               
                                                                               
        self.Graph = gclass.Graph(self.graph)                                  
                                                                               
        not_alone_atoms = []                                                   
                                                                               
        for key in self.graph:                                                 
            degree = self.Graph.vertex_degree(key)                             
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
                    path = self.Graph.find_all_paths(key, key2)                
                    if len(path) == 1:                                         
                        if len(path [0]) == 2:                                 
                            
                            print(key+1, key2+1, len(path),len(path [0]), self.bond_type[(key,key2)]  )     #
                            self.rotatable_bonds.append([key,key2])
        print self.rotatable_bonds
        
        for rotatable_bond in self.rotatable_bonds:
            self.Graph.find_rotatable_blocks(rotatable_bond[0], rotatable_bond[1])
        #-----------------------------------------------------------------------------------------------------------#
        
        
        
        
        


        
        
        
         




class Atom:
    """ Class doc """
    
    def __init__ (self, name         ='Xx',
                        index        =None, 
                        symbol       ='X', 
                        pos          = None, 
                        resi         = None, 
                        resn         = None, 
                        chain        = '', 
                        atom_id      = 0, 
                        molecule     = None,
                        #Vobject_id   = None, 
                        #Vobject_name = '', 
                        #Vobject      = None
                        ):
 
        if pos is None:
            pos = np.array([0.0, 0.0, 0.0])

        self.pos      = pos     # - coordinates of the first frame
        self.index    = index   # 
        self.name     = name    #
        self.symbol   = symbol  #
        self.resi     = resi    #
        self.resn     = resn    #
        self.chain    = chain   #
        self.molecule = molecule
        #self.Vobject = Vobject
        #self.residue = residue    

        #self.atom_id  = atom_id        # An unique number
        #self.color    = at.get_color    (name)
        #self.color_id = None
        #self.color.append(0.0)  

        #self.col_rgb      = at.get_color_rgb(name)
        #self.radius       = at.get_radius   (name)
        #self.vdw_rad      = at.get_vdw_rad  (name)
        #self.cov_rad      = at.get_cov_rad  (name)
        #self.ball_radius  = at.get_ball_rad (name)
        self.connected      = []
