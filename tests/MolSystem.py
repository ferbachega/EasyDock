
class Molecule:
    """ Class doc """
    
    def __init__ (self, atoms    = []   , 
                        bonds    = []   , 
                        name     = 'UNK', 
                        index    = None ,
                        #chain   = None,
                        #Vobject = None
                        ):
        """ Class initialiser """
        self.atoms      = atoms
        self.bonds      = bonds
        
        self.bond_types = {}
        
        self.name       = name
        self.index      = index
        
        
        #self.resi     = index
        #self.resn     = name
        #self.chain    = chain
        #self.Vobject  = Vobject



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
