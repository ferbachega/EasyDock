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
        
        # building the graph 
        self.graph = {}                                                        
        for atom in self.atoms:                                                
            self.graph[atom.index] = []                                        
                                                                               
        for bond in self.index_bonds_pairs:                                    
            self.graph [ bond[0]].append(bond[1])                              
            self.graph [ bond[1]].append(bond[0])                              
                                                                      
        self.Graph = gclass.Graph(self.graph)                                  

        #------------------------------------------------------------------------------------------------------------#
        # finding atoms with only one connection                                                                        
        not_alone_atoms = []                                                   
                                                                               
        for key in self.graph:                                                 
            degree = self.Graph.vertex_degree(key)                             
            print(key, degree)                                                 
            if degree > 1:                                                     
                not_alone_atoms.append(key)                                    
        #------------------------------------------------------------------------------------------------------------#
        
        
        #------------------------------------------------------------------------------------------------------------#
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
                            
                            #for  atom in 
                            
                            #print(key+1, key2+1, len(path),len(path [0]), self.bond_type[(key,key2)]  )     #
                            self.rotatable_bonds.append([key,key2])
        print self.rotatable_bonds
        
        
        
        
        
        # --------------------------------------------------------
        '''
        finding the methyl groups
        '''
        _buffer = []
        
        methyl_groups = []
        for bond in self.rotatable_bonds:
            for vertex_v in bond:
                if vertex_v in _buffer:
                    pass
                else:
                    _buffer.append(vertex_v)
                    n = 0 
                    
                    for vertex_w in self.graph[vertex_v]:
                        
                        if self.atoms[vertex_w].name == 'H':
                            print vertex_v, vertex_w, 'H'
                            n += 1
                    
                    if n == 3:
                        print vertex_v, 'methyl'
                        methyl_groups.append(vertex_v)
        # --------------------------------------------------------
        
        
        
        # --------------------------------------------------------
        '''excluding methyl groups from rotatable bonds'''

        self.rotatable_bonds_edited  = []
        for bond in self.rotatable_bonds:
            if bond[0] in methyl_groups or bond[1] in methyl_groups:
                pass
            else:
                self.rotatable_bonds_edited.append(bond)
        # --------------------------------------------------------
               
        
        
        
        
        
        
        #for rotatable_bond in self.rotatable_bonds_edited:
        #    self.Graph.find_rotatable_blocks(rotatable_bond[0], rotatable_bond[1])
        #-----------------------------------------------------------------------------------------------------------#
        
        #_ATOMLINEFORMAT1    = "{:<6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4s}{:2s}  {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}{:2s}{:2s}"

        self.Graph.find_rotatable_blocks2(rot_bonds = self.rotatable_bonds_edited)
        
        self.block_connections = {}
        for block in self.Graph.blocks:
            self.block_connections[block] = []
        
        
        for bond in self.rotatable_bonds_edited:
            
            blocka = None
            blockb = None 
            
            for block in self.Graph.blocks:
                if bond[0] in self.Graph.blocks[block]:
                    #print bond[0], block,  self.Graph.blocks[block]
                    blocka = block
                if bond[1] in self.Graph.blocks[block]:
                    #print bond[1], block,  self.Graph.blocks[block]
                    blockb = block
                #print bond, block, self.Graph.blocks[block]
            
            #print  bond[0], bond[1], blocka, blockb
            self.block_connections[blocka].append([blockb,[bond[0], bond[1]]])
            self.block_connections[blockb].append([blocka,[bond[1], bond[0]]])
            #self.block_connections[blocka] = [blockb, [bond[0], bond[1]]]
            #self.block_connections[blockb] = [blocka, [bond[1], bond[0]]]
        
        print self.block_connections
        #self.export_PDB (fileout = 'molecule.pdb')
        
        
    def export_PDB (self, fileout = 'molecule.pdb'):
        """ Function doc """
        fileout = open(fileout, 'w')
        
        lines = []
        _ATOMLINEFORMAT1    = "{:<6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4s}{:2s}  {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}{:2s}{:2s}\n"
        for block in self.Graph.blocks:
            lines.append('REMARK Block '+ str(block)  +' \n')
            print 'REMARK Block ', block, '\n'
            
            for index in self.Graph.blocks[block]:
                at = self.atoms[index]
                #print "ATOM  ", index, self.atoms[index].name,  self.atoms[index].pos
                lines.append(_ATOMLINEFORMAT1.format ( "ATOM", index + 1, at.name , " ", 'UNK', 'A' , '', str(block),
                                                  
                                                  at.pos[0], at.pos[1], at.pos[2], 0.00, at.charge, '', at.name, " " ))
                print ( _ATOMLINEFORMAT1.format ( "ATOM", index + 1, at.name , " ", 'UNK', 'A' , '', str(block),
                                                  
                                                  at.pos[0], at.pos[1], at.pos[2], 0.00, at.charge, '', at.name, " " ) )
        
        fileout.writelines(lines)
        


class Block:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        self.atoms = []
        
        self.connections = []
        
        self.bonds
        
        pass


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
                        charge       = 0.00,
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
        self.charge   = charge
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
