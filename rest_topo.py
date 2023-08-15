#!/bin/env python
import numpy as np
import pandas as pd
​
class load_topology(filename):
   def __init__(self):
      self.input_at = None
      scaled_sections = li st
      molecules = dict
      sections = dict
      hard_order_sections = [ "defaults", "atomtypes", "nonbond_params", "bondtypes", "constrainttypes", ]
      with open(filename) as topo:
         self.readfile = topo.readlines():
​
   def parse_section(self,trunk_list):
      first_round = True
      output = []
      for line in trunk_list:
         if "[" not in line:
            output.append(line)
         if ";" == line[0]:
            #output.append(line)
            continue
         if "[" in line: 
            break
      return output
​
   def parse_molecules(self):
      is_molecule = [ i for i, line in enumerate(self.readfile) if "moleculetype" in line]
      
​
   def gather_atomtypes(self):
      is_atomtypes = [ i for i, line in enumerate(self.readfile) if "atomtypes" in line]
      in_atomtypes = [" [ atomtypes ] ", "; name      at.num  mass     charge ptype  sigma      epsilon"]
      for i in range(len(is_atomtypes)-1):
         in_atomtypes += parse_section[is_atomtypes[i]:is_atomtypes[i+1]]
      in_atomtypes += parse_section[is_atomtypes[-1]:]