#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  MOL2Files.py
#  
#  Copyright 2018 Fernando Bachega <fernando@Fenrir>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
from pprint import pprint
import numpy as np
import MolSystem
import MolecularSystem


'''
@<TRIPOS>MOLECULE
DCM Pose 1
   32    33     0     0     0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 N          63.2680   27.8610   32.2290 N.3     4  VAL4        0.0000

      1 C1         18.8934    5.5819   24.1747 C.2       1 <0>       -0.1356 
      2 C2         18.1301    4.7642   24.8969 C.2       1 <0>       -0.0410 
      3 C3         18.2645    6.8544   23.7342 C.2       1 <0>        0.4856 
      4 C4         16.2520    6.2866   24.7933 C.2       1 <0>        0.8410 
      5 C5         15.3820    3.0682   25.1622 C.3       1 <0>        0.0000 
      6 C6         15.4162    1.8505   26.0566 C.3       1 <0>        0.2800 
      7 C7         16.7283    2.0138   26.8111 C.3       1 <0>        0.2800 
      8 C8         16.0764    4.1199   26.0119 C.3       1 <0>        0.5801 
      9 C9         17.9106    1.3823   26.0876 C.3       1 <0>        0.2800 
     10 N1         17.0289    7.1510   24.0411 N.2       1 <0>       -0.6610 
     11 N2         16.8196    5.0644   25.2302 N.am      1 <0>       -0.4691 
     12 N3         19.0194    7.7275   22.9859 N.pl3     1 <0>       -0.8500 
     13 O1         18.7676   -2.3524   26.1510 O.3       1 <0>       -1.0333 
     14 O2         20.3972   -0.3812   26.2318 O.3       1 <0>       -1.0333 
     15 O3         15.0888    6.5824   25.0727 O.2       1 <0>       -0.5700 
     16 O4         18.9314   -0.7527   24.1606 O.2       1 <0>       -1.0333 
     17 O5         16.9690    3.4315   26.8994 O.3       1 <0>       -0.5600 
     18 O6         14.3223    1.8946   26.9702 O.3       1 <0>       -0.6800 
     19 O7         17.9091   -0.0135   26.3390 O.3       1 <0>       -0.5512 
     20 P1         19.0969   -0.9440   25.6653 P.3       1 <0>        1.3712 
     21 H1         19.9176    5.3550   23.9105 H         1 <0>        0.1500 
     22 H2         18.5100    3.8155   25.2595 H         1 <0>        0.1500 
     23 H3         15.8520    2.8983   24.1870 H         1 <0>        0.0000 
     24 H4         14.3405    3.3601   24.9711 H         1 <0>        0.0000 
     25 H5         15.3663    0.9351   25.4839 H         1 <0>        0.0000 
     26 H6         16.6681    1.6130   27.8171 H         1 <0>        0.0000 
     27 H7         15.3483    4.6961   26.6094 H         1 <0>        0.0000 
     28 H8         18.8490    1.8078   26.4511 H         1 <0>        0.0000 
     29 H9         17.8303    1.5497   25.0110 H         1 <0>        0.0000 
     30 H10        19.9527    7.4708   22.7715 H         1 <0>        0.4000 
     31 H11        18.5977    8.5756   22.6932 H         1 <0>        0.4000 
     32 H12        14.2530    1.0535   27.4278 H         1 <0>        0.4000 
@<TRIPOS>BOND
    1     1     2 2
    2     1     3 1
    3     2    11 1
    4     3    10 2
    5     3    12 1
    6     4    10 1
    7     4    11 am
    8     4    15 2
    9     5     6 1
   10     5     8 1
   11     6     7 1
   12     6    18 1
   13     7     9 1
   14     7    17 1
   15     8    11 1
   16     8    17 1
   17     9    19 1
   18    13    20 1
   19    14    20 1
   20    16    20 2
   21    19    20 1
   22     1    21 1
   23     2    22 1
   24     5    23 1
   25     5    24 1
   26     6    25 1
   27     7    26 1
   28     8    27 1
   29     9    28 1
   30     9    29 1
   31    12    30 1
   32    12    31 1
   33    18    32 1
@<TRIPOS>SUBSTRUCTURE
'''





def get_atom_list_from_mol2_frame (raw_atoms, frame = True, gridsize = 3):
    """ Function doc """
    #nCPUs =  multiprocessing.cpu_count()
    #pool  = multiprocessing.Pool(nCPUs)
    #pdb_file_lines  = frame.split('\n')   
    #atoms = (pool.map(parse_pdb_line, pdb_file_lines))

    atoms             = []
    frames            = []
    frame_coordinates = []
    
    for line in raw_atoms:
    
        line = line.split()
        index  = 0
        if len(line) > 1:
            #print (line) 
            index    = int(line[0])-1
            
            at_name  = line[1]
            
            at_pos   = np.array([float(line[2]), float(line[3]), float(line[4])])
            
            at_resi  = int(line[6])
                     
            at_resn  = line[6]
                     
            at_ch    = 'X'          
            

            at_symbol = line[5].split('.')
            at_symbol = at_symbol[0]

            at_charge = float(line[-1])*100

            gridpos   = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
            
            #atoms.append([index, at_name,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos ])
            
            #atom =  MolSystem.Atom (name         = at_name,
            #                        index        = index, 
            #                        symbol       = at_symbol, 
            #                        pos          = at_pos, 
            #                        resi         = at_resi, 
            #                        resn         = at_resn, 
            #                        chain        = at_ch, 
            #                        atom_id      = 0, 
            #                        molecule     = None,
            #                        charge       = at_charge
            #                        #Vobject_id   = None, 
            #                        #Vobject_name = '', 
            #                        #Vobject      = None
            #                        )
            
            
            atom =  MolecularSystem.Atom(index        = index    , 
                                         name         = at_name  ,
                                         symbol       = at_symbol, 
                                         atom_type    = None     , 
                                         pos          = at_pos   , 
                                         charge       = at_ch    ,
                                         
                                         sector       = None,
                                         resi         = None, 
                                         resn         = None, 
                                         chain        = None, 
                                         
                                         bonds        = []  ,
                                         molecule     = None,
                                         )
            
            

            atoms.append(atom)
            index += 1
    #print (atoms)
    return atoms

def get_bonds (raw_bonds, atoms):
    """ Function doc """
    index_bonds              = []
    index_bonds_pairs        = []
    index_bonds_pairs_orders = []
    
    bonds = {}
    
    bond_type = {}
    
    
    
    #print (raw_bonds)
    print ('Obtain bonds from original MOL2 file')
    for line in raw_bonds:
        line = line.split()
        #print line 
        
        if line == []:
            pass
        else:
            if len(line) == 4:
                #print line
                index    = int(line[0])-1           
                atom1    = int(line[1])-1
                atom2    = int(line[2])-1
                order    = line[3]
                
                index_bonds      .append(atom1)
                index_bonds      .append(atom2)
                
                bonds[index] = [atom1, atom2, order]
                index_bonds_pairs.append([atom1,atom2, order])
                
                bond_type[(atom1, atom2)] = order
                bond_type[(atom2, atom1)] = order
                
        '''
        if len(line) == 4:
            index    = int(line[0])            
            atom1    = int(line[1]-1)
            atom2    = int(line[2]-1)
            order    = line[3]

            index_bonds      .append(atom1)
            index_bonds      .append(atom2)
            index_bonds_pairs.append([atom1,atom2])
            
            index_bonds_pairs_orders.append(order)
        '''
    #print bond_type
    #print index_bonds_pairs
    #return [bonds, index_bonds, index_bonds_pairs]
    
    Bonds = []
    n = 0
    for raw_bond in index_bonds_pairs:
        index_a = raw_bond[0]
        index_b = raw_bond[1]
        b_type  = raw_bond[2]
        bond =  MolecularSystem.Bond(index    = n             , 
                                     atom1    = atoms[index_a], 
                                     atom2    = atoms[index_b], 
                                     b_type   = b_type        , 
                                     sector   = None, 
                                     molecule = None)
        Bonds.append(bond)
        n += 1
    
    return Bonds #bonds, index_bonds, index_bonds_pairs, bond_type











def parser_raw_mol2_info (raw_molecule, log = False):
    """ Function doc """
    firstmolecule =  raw_molecule.split('@<TRIPOS>ATOM')
    header        =  firstmolecule[0]
    firstmolecule =  firstmolecule[1].split('@<TRIPOS>BOND')
    raw_atoms     =  firstmolecule[0]
    raw_bonds     =  firstmolecule[1]
   
    header    = header.split('\n')
    #print header
    mol_name   = header[1]
    mol_size   = header[2]
    mol_type   = header[3]
    mol_charge = header[4]
    mol_info   = header[5]
    
    raw_atoms = raw_atoms.split('\n')
    raw_bonds = raw_bonds.split('\n')
    
    atoms = get_atom_list_from_mol2_frame(raw_atoms = raw_atoms, frame = True,  gridsize = 3)

    #bonds, index_bonds, index_bonds_pairs, bond_type = get_bonds(raw_bonds, atoms)
    Bonds = get_bonds(raw_bonds, atoms)

    
    #print 'Bonds', Bonds
    #print 'Atoms', atoms
    #print '\n\n'
    #print index_bonds
    #print '\n\n'
    #print index_bonds_pairs
    #print '\n\n'

    #molecule = MolSystem.Molecule ( name              = mol_name  ,
    #                                size              = mol_size  ,
    #                                _type             = mol_type  ,
    #                                charge            = mol_charge,
    #                                info              = mol_info  ,
    #                                index             = None         ,
    #
    #                                atoms             = atoms     , 
    #                                bonds             = bonds     , 
    #                                index_bonds_pairs = index_bonds_pairs ,
    #                                bond_type         = bond_type         ,
    #                                )
    #
    
    print 'Building molecule'
    Molecule =  MolecularSystem.Molecule(index       = None     , 
                                         name        = mol_name ,
                                         atoms       = atoms    , 
                                         bonds       = Bonds    , 
                                         sectors     = []       ,
                                         connections = []       
                                         )
    

    Molecule.building_graph()
    #Molecule.printAtoms()
    #Molecule.printBonds()
    print 'Molecule Done'
    
    
    
    return Molecule


def load_mol2_files (infile = None, VMSession =  None, gridsize = 3):
    """ Function doc """
    print ('\nstarting: parse_mol2')

    #initial = time.time()
    
    
    molecule_list = []
    
    
    with open(infile, 'r') as mol2_file:
        pdbtext = mol2_file.read()

        raw_molecules     =  pdbtext.split('@<TRIPOS>MOLECULE')
        
        number_of_molecules = len(raw_molecules)
        print ('Number of molecules:', number_of_molecules)
        
        n = 1
        for raw_molecule in raw_molecules:
            #print ('molecule', n)
            #print raw_molecule
            
            if '@<TRIPOS>ATOM' in raw_molecule:
                
                
                
                print ('molecule', n, True )
                molecule = parser_raw_mol2_info(raw_molecule, log = False)
                molecule_list.append(molecule)
            
            else: 
                print ('molecule', n, False)
            n += 1            
    
               
    return molecule_list


def load_file (infile = None, EasyDockSession =  None, gridsize = 3):
    """ Function doc """
    print ('\nstarting: parse_mol2')
    molecule_list = []
    
    
    
    with open(infile, 'r') as mol2_file:
        pdbtext = mol2_file.read()

        raw_molecules       =  pdbtext.split('@<TRIPOS>MOLECULE')
        number_of_molecules = len(raw_molecules)
        print ('Number of molecules:', number_of_molecules)
        n = 1
        
        for raw_molecule in raw_molecules:        
            if '@<TRIPOS>ATOM' in raw_molecule:
                
                print ('molecule', n, True )
                molecule = parser_raw_mol2_info(raw_molecule, log = False)
                molecule_list.append(molecule)
            
            else: 
                print ('molecule', n, False)
            n += 1            
    
               
    return molecule_list







#load_mol2_files (infile = '/home/fernando/programs/EasyDock/mol2/com7.mol2', VMSession =  None, gridsize = 3)

