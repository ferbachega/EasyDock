#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ECore.py
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

import MOL2Files_Reader as Mol2Files


class EasyDock:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
        
        self.ligands  = []
        self.receptor = []
        
        
    
    def LoadLigandToSystem (self, filein = None, _type = 'mol2'):
        """ Function doc """
        if filein == None:
            pass
        
        else:
        
            if _type ==  'mol2':
                self.LoadMol2FileToSystem (filein, _type = 'ligand')
        


    def LoadMol2FileToSystem (self, filein, _type = 'ligand'):
        """ Function doc """
        molecules = Mol2Files.load_file(filein)
        
        if _type == 'ligand':
            for molecule in molecules:
                self.ligands.append(molecule)

        else:
            for molecule in molecules:
                self.receptor.append(molecule)
        
        


easydock =  EasyDock()

easydock.LoadLigandToSystem(filein = '/home/fernando/programs/EasyDock/mol2/com7.mol2')
        
        
