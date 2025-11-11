#!/usr/bin/env python3



from numpy import exp, log, arange, array, float32, unique, sqrt, where, ndarray, isin, dtype
from copy import deepcopy


def generate_kappas(nreps:int, kappa_max):
   return   exp( arange(nreps) * log(kappa_max) / (nreps-1))
      
def generate_temperature_ladder(nreps:int, tlow, thigh):
   return tlow * exp( arange(nreps) * log(thigh/tlow) / (nreps-1))

def compute_lambda(temperatures:array):
   return temperatures.min()/temperatures

def compute_se(e1:float,e2:float,s1:float,s2:float):
   eij = (float32(e1)*float32(e2))**0.5
   sij = 0.5*(float32(s1)+float32(s2))
   return eij,sij

def generate_nonbonded(at1, at2):
   e1 = float32(at1[-1])
   e2 = float32(at2[-1])
   s1 = float32(at1[-2])
   s2 = float32(at2[-2])
   funct = str(1)
   eij, sij = compute_se(e1,e2,s1,s2)
   a1, a2 = at1[0], at2[0]
   return f'{a1:>5} {a2:>4} {funct:^5} {sij:>10.4f} {eij:>8.4f}\n'

def query_typed(query, qtype):
   return eval(qtype)(input(query))

def molecule_selections(queryN:list, queryM:list):
   response = query_typed(*queryN)
   molecule_sel = []
   for i in range(response):
      molecule_sel.append(query_typed(*queryM))
   return molecule_sel

def _parse_section(readfile,trunks):
   first_round = True
   output = []
   for line in readfile[trunks:]:
      if "[" not in line and first_round!=True and len(line.split()) != 0:
         output.append(line)
      elif ";" == line[0]: 
         continue
      elif "[" in line and first_round!=True: 
         break
      first_round = False
   return output

def get_molecule_atomtypes(sections,molecule:int=0):
   atomtype_list = []
   for i in sections['moleculetype'][molecule]['atoms']:
      if len(i.split())>1 and ';' not in i.split()[0] and '[' not in i :
         atomtype_list.append(i.split()[:2])
   dataset = array(atomtype_list,dtype=object)
   atomtypes = dataset[:,1]
   return unique(atomtypes)

def dihedraltype_array(dihedraltypes:list):
   data = []
   dt = dtype([('i','U10'),
                  ('j', 'U10'),  # 2nd string
                  ('k', 'U10'),  # 3rd string
                  ('l', 'U10'),  # 4th string
                  ('func', 'i4'),    # 5th integer
                  ('angle', 'f4'),    # 6th float
                  ('K', 'f4'),    # 7th float
                  ('mult', 'i4')])   # 8th integer])
   for string in dihedraltypes:
      line = string[:string.find(';')].split()
      data.append((line[0], line[1], line[2], line[3], int(line[4]), float(line[5]), float(line[6]), int(line[7])))
   return array(data, dtype=dt)

def dihedral_array(dihedraltypes:ndarray, dihedrals:list, atom_mapping:dict, verbose=False):
   data = []
   search_dihedrals = []
   dtype_dihedrals = dtype([('i','i4'),
                  ('j', 'i4'),  # 2nd string
                  ('k', 'i4'),  # 3rd string
                  ('l', 'i4'),  # 4th string
                  ('func', 'i4'),    # 5th integer
                  ('angle', 'f4'),    # 6th float
                  ('K', 'f4'),    # 7th float
                  ('mult', 'i4')])   # 8th integer])
   for string in dihedrals:
      line = string[:string.find(';')].split()
      if len(line) == 5:
         search_dihedrals.append(tuple(map(int,line)))
      elif len(line) == 8:
         data.append((int(line[0]), int(line[1]), int(line[2]), int(line[3]), int(line[4]), float(line[5]), float(line[6]), int(line[7])))
      elif not line:
         continue
      else:
         raise RuntimeError(f'Line not processed: {line}')
   print(search_dihedrals)
   for (i, j, k, l, func) in search_dihedrals:
      if func == 4:
         mask = ((isin(dihedraltypes['i'], atom_mapping[i]) & \
               isin(dihedraltypes['j'], atom_mapping[j]) & \
               isin(dihedraltypes['k'], atom_mapping[k]) & \
               isin(dihedraltypes['l'], atom_mapping[l]) & \
               isin(dihedraltypes['func'], func)) | \
               (isin(dihedraltypes['i'], atom_mapping[l]) & \
               isin(dihedraltypes['j'], atom_mapping[k]) & \
               isin(dihedraltypes['k'], atom_mapping[j]) & \
               isin(dihedraltypes['l'], atom_mapping[i]) & \
               isin(dihedraltypes['func'], func))) & isin(dihedraltypes['func'], func)
         mask_Xl = (isin(dihedraltypes['i'], 'X') & \
               isin(dihedraltypes['j'], atom_mapping[j]) & \
               isin(dihedraltypes['k'], atom_mapping[k]) & \
               isin(dihedraltypes['l'], atom_mapping[l]) & \
               isin(dihedraltypes['func'], func)) 
         mask_Xll = (isin(dihedraltypes['i'], 'X') & \
               isin(dihedraltypes['j'], 'X') & \
               isin(dihedraltypes['k'], atom_mapping[k]) & \
               isin(dihedraltypes['l'], atom_mapping[l]) & \
               isin(dihedraltypes['func'], func))
         if verbose:
            if any(mask):
               print(f'found complete match for {(i,j,k,l,func)}')
               print(f'{" ".join([atom_mapping[atom] for atom in [i, j, k, l]])}')
               print(dihedraltypes[mask])
            elif any(mask_Xl):
               print(f'found match three for {(i, j, k, l,func)}')
               print(f'{" ".join([atom_mapping[atom] for atom in [i, j, k, l]])}')
               print(dihedraltypes[mask_Xl])
            elif any(mask_Xll):
               print(f'found match two for {(i, j, k, l,func)}')
               print(f'{" ".join([atom_mapping[atom] for atom in [i, j, k, l]])}')
               print(dihedraltypes[mask_Xll])
            else:
               raise LookupError(f'Dihedral not found: {[i, j, k, l, func].__str__()} {[atom_mapping[atom] for atom in [i, j, k, l]].__str__()}')
         elif not any(mask | mask_Xl | mask_Xll):
            raise LookupError(f'Dihedral not found: {[i, j, k, l, func].__str__()} {[atom_mapping[atom] for atom in [i, j, k, l]].__str__()}')
         
         if any(mask):
            for entries in dihedraltypes[mask]:
               data.append((i, j, k, l, func, entries[5], entries[6], entries[7])) 
         elif any(mask_Xl):
            for entries in dihedraltypes[mask_Xl]:
               data.append((i, j, k, l, func, entries[5], entries[6], entries[7])) 
         else:
            for entries in dihedraltypes[mask_Xll]:
               data.append((i, j, k, l, func, entries[5], entries[6], entries[7]))
            
      if func == 9:
         mask = ((isin(dihedraltypes['i'], atom_mapping[i]) & \
               isin(dihedraltypes['j'], atom_mapping[j]) & \
               isin(dihedraltypes['k'], atom_mapping[k]) & \
               isin(dihedraltypes['l'], atom_mapping[l]) & \
               isin(dihedraltypes['func'], func)) | \
               (isin(dihedraltypes['i'], atom_mapping[l]) & \
               isin(dihedraltypes['j'], atom_mapping[k]) & \
               isin(dihedraltypes['k'], atom_mapping[j]) & \
               isin(dihedraltypes['l'], atom_mapping[i]) & \
               isin(dihedraltypes['func'], func))) & isin(dihedraltypes['func'], func)
         mask_Xlr = ((isin(dihedraltypes['i'], 'X') & \
               isin(dihedraltypes['j'], atom_mapping[j]) & \
               isin(dihedraltypes['k'], atom_mapping[k]) & \
               isin(dihedraltypes['l'], 'X') & \
               isin(dihedraltypes['func'], func)) | \
               (isin(dihedraltypes['i'], 'X') & \
               isin(dihedraltypes['j'], atom_mapping[k]) & \
               isin(dihedraltypes['k'], atom_mapping[j]) & \
               isin(dihedraltypes['l'], 'X') & \
               isin(dihedraltypes['func'], func))) & isin(dihedraltypes['func'], func)
         if verbose:
            if any(mask):
               print(f'found complete match for {(i,j,k,l,func)}')
               print(f'{" ".join([atom_mapping[atom] for atom in [i, j, k, l]])}')
               print(dihedraltypes[mask])
            elif any(mask_Xlr):
               print(f'found match two for {[i, j, k, l]}')
               print(f'{" ".join([atom_mapping[i] for i in [i, j, k, l]])}')
               print(dihedraltypes[mask_Xlr])
            else:
               raise LookupError(f'Dihedral not found: {[i, j, k, l, func].__str__()} {[atom_mapping[atom] for atom in [i, j, k, l]].__str__()}')
         elif not any(mask | mask_Xlr):
            raise LookupError(f'Dihedral not found: {[i, j, k, l, func].__str__()} {[atom_mapping[atom] for atom in [i, j, k, l]].__str__()}')
         if any(mask): 
            for entries in dihedraltypes[mask]:
               data.append((i, j, k, l, func, entries[5], entries[6], entries[7])) 
         else:
            for entries in dihedraltypes[ mask_Xlr]:
               data.append((i, j, k, l, func, entries[5], entries[6], entries[7])) 
      
   
   structured_array = array(data, dtype=dtype_dihedrals)
   structured_array['func'] *= -1
   structured_array.sort(order=['l','k','j','mult','i'])
   structured_array.sort(order='func')
   structured_array['func'] *= -1
   
   return structured_array

def dihedraltypes_strings_list(dihedraltypes):
   stringsout = [f'{dihedraltype[0]:<5} {dihedraltype[1]:<5} {"s"+dihedraltype[2]:<5} {"s"+dihedraltype[3]:<5} {dihedraltype[4]:^9}' + \
                 f'{dihedraltype[5]:<10}{dihedraltype[5]:<10.5f}{dihedraltype[7]}\n' \
                 for dihedraltype in dihedraltypes]
   return stringsout

def dihedrals_strings_list(dihedrals):
   stringsout = [f'{dihedral[0]:<5} {dihedral[1]:<5} {dihedral[2]:<5} {dihedral[3]:<5} {dihedral[4]:^9}' + \
                 f'{dihedral[5]:<10.3f}{dihedral[6]:<10.5f}{dihedral[7]}\n' \
                 if len(dihedral) == 8 else f'{dihedral[0]:<5} {dihedral[1]:<5} {dihedral[2]:<5} {dihedral[3]:<5} {dihedral[4]:^9}' \
                 for dihedral in dihedrals ]
   return stringsout
 
def topology_writer(**kwargs):
   ofile = kwargs.get('outfile', 'topol')
   filepath = kwargs.get('filepath', './')
   lambda_i = kwargs.get('lambda_i', None)
   kappa_i = kwargs.get('kappa_i', None)
   lines = kwargs.get('lines', None)
   verbose = kwargs.get('verbose', False)
   if any(lines):
      fullpath = f'{filepath}{ofile}'
      for additional in [lambda_i, kappa_i]: fullpath = f'{fullpath}-{additional:0.3f}' if additional else fullpath
      fullpath += ".top"
      with open(fullpath,'w') as file:
         file.writelines(lines)
      if verbose: print(f'Saved {fullpath}')
   else: print('try again with self.run()')
   

class topo2rest():
   def __init__(self, ifile:str='processed.top', temps:list = [300.0, 500.0], nreps:int = 20, **kwargs):
      '''Convert processed topology to REST2/3 input topology
         E_tot = gamma*E^{pp} + sqrt(gamma)*E^{pw} + E^{ww}
         gamma = T_0/T_i
         for REST2/3 bonds and angles are not scaled per the REST2 paper,
         however, for REST2 LJ parameters epsilon_i is scaled by epsilon_i*gamma
         for tempered atoms and all others are unmodified. 
         In the case of REST3, with the additon of sqrt(gamma)*kappa*E^{pw} which differs 
         from sqrt(gamma)*E^{pw}, we produced the combination rule 2 nonbonded terms
         involving protein-water interactions to override the pw nonbonded interactions
         as including gamma*kappa with each hot atom epsilon_i would result in the 
         incorrect form of: E_tot = gamma*kappa*E^{pp} + sqrt(gamma*kappa)*E^{pw} + E^{ww}:
         ra`ther than the correct form: 
         E_tot = gamma*kappa*E^{pp} + sqrt(gamma)*kappa*E^{pw} + E^{ww}
         input
         ifile = inputtopology.top
         temps = ["lower temp":float, "upper temp":float ]; temperature range of replicas
         kappa:float = kappa scaling; if not equal to 1 REST3 implementation active
         hot_molecules = hot molecules; [0] for most systems will select the protein'''

      self.nreps = nreps
      self._processed = False
      self._lambdai = None
      self._kappa = None
      self.kappa_atoms = None
      self.hot_molecules = None
      self.nmol = None
      self.molecules = None
      self.templadder = generate_temperature_ladder(nreps, *temps)
      self.lambdai = self.templadder
      self._sections = {}
      self._sections_out = {}
      self._scaled_dihedrals = {}
      self._scaled_atomtypes = {}
      self._scaled_molecules={}
      self._scaled_nonbonded={}
      self._molecule_atoms = {}
      self.kappa_low_temp = 330
      self.kappa = 1.0
      self._methods = { 'ssrest3':['lambda','kappa'], 'rest2':['lambda'], 'solvent_scaled':['kappa'] }
      
      self.hard_order_sections = [ "defaults", "atomtypes", "nonbond_params", "bondtypes", \
                                   "constrainttypes", "angletypes", "dihedraltypes","system", "molecules", "moleculetype"]
      with open(ifile) as topo:
         self.readfile = topo.readlines()
      self._gather_param_sections()
      self._atomtypes = [i[:i.find(';')].split() \
                               for i in self._sections['atomtypes'] if ';' not in i[:3] if '[' not in i[:3] \
                               if '\n' not in i[:3]]
      self._ordered_molecules()
      self.molecule_atoms = True
      if kwargs.get('verbose', False):
         print('topology file read')
   
   @property
   def lambdai(self):
      return self._lambdai
   
   @lambdai.setter
   def lambdai(self, temps):
      if isinstance(temps, int):
         self._lambda = arange(1,temps+1)
      else:
         self._lambdai = compute_lambda(temps)
      pass
       
   @property
   def kappa(self):
      return self._kappa
   
   @kappa.setter
   def kappa(self, kappa_max):
      assert isinstance(kappa_max, (int, float)), 'kappa_max must be integer or float'
      kappa = generate_kappas(self.nreps, kappa_max=kappa_max)
      kappa[where(self.templadder<=self.kappa_low_temp)[0]] = 1.0
      self._kappa = {i:j for i,j in \
                     zip(self.lambdai,kappa)}
      pass
      
   @property
   def molecule_atoms(self):
      return self._molecule_atoms
   
   @molecule_atoms.setter
   def molecule_atoms(self, a):
      if a == True:
         atdict = {}
         for nmol in range(len(self._sections['moleculetype'])):
            atdict[nmol] = {int(i.split()[0]): i.split()[1] for i in self._sections['moleculetype'][nmol]['atoms'] if len(i.split()) != 0 and ';' not in i[:5]}
         self._molecule_atoms = atdict
      pass
      
   def _get_scaled_molecules(self, lambda_on=True):
      for lambdai in self.lambdai:
         self._scaled_molecules[lambdai]=deepcopy(self._sections['moleculetype'])
         if lambda_on:
            for hot in self.hot_molecules:
               atoms=[i.split() if len(i.split()) == 8 else i.split()[:-3] for i in self._sections['moleculetype'][hot]['atoms'] if len(i.split()) != 0]
               scaled_charges=[]
               for i in atoms:
                  s_charge=float(i[6])*sqrt(lambdai) # scalling with sqrt(lambda)
                  stringout = f'{i[0]:>6} {"s"+i[1]:>10} {i[2]:>6} {i[3]:>6} {i[4]:>6} {i[5]:>6} {s_charge:>10.6e} {float(i[7]):>10.3f}\n'
                  scaled_charges.append(stringout)
               self._scaled_molecules[lambdai][hot]['atoms'] = scaled_charges
               self._scaled_molecules[lambdai][hot]['dihedrals'] = deepcopy(self._scaled_dihedrals[hot][lambdai])
            
   def _get_scale_dehedrals(self, lambda_on=True):
      
      dihedral_types = [ i for i in self._sections['dihedraltypes'] if len(i[:i.find(';')].split()) > 4 ]
      if lambda_on:
         for hot in self.hot_molecules:
            dihedrals = dihedral_array(dihedraltype_array(dihedral_types),self._sections['moleculetype'][hot]['dihedrals'],\
                                                self.molecule_atoms[0])
            unscaled_K = deepcopy(dihedrals['K'])
            self._scaled_dihedrals[hot] = {}
            for lambdai in self.lambdai:
               dihedrals['K'] = unscaled_K * float32(lambdai)
               self._scaled_dihedrals[hot][lambdai] = dihedrals_strings_list(dihedrals)      
      else:
         for hot in self.hot_molecules:
               self._scaled_dihedrals[hot] = {i:self._sections['moleculetype'][hot]['dihedrals'] for i in self.lambdai}
      pass

   def _populate_out(self, **kwargs):
      if kwargs.get('verbose', False):
         method = kwargs.get('method', None)
         if not method == 'solvent_scaled':
            print(f'Temperature Ladder:         {" ".join(self.templadder)}')
            print(f'Corresonding Lambdas:       {" ".join(self.lambdai)}')
         if kwargs.get('method',None) in [i for i in self._methods.keys() if i != 'rest2']:
            print(f'Corresponding kappa values: {" ".join(self.kappa)}')
         print(f'Produce {self.nreps} replicas')
      for lambda_i in self.lambdai:
         lines=self._sections['defaults']+['\n']
         lines+=['[ atomtypes ]\n']+self._scaled_atomtypes[lambda_i]+['\n']
         lines+=['[ nonbond_params ]\n']+self._scaled_nonbonded[lambda_i]+['\n']
         lines+=self._sections['bondtypes']+['\n']
         lines+=self._sections['constrainttypes']+['\n']
         lines+=self._sections['angletypes']+['\n']
         lines+=self._sections['dihedraltypes']+['\n']
         for molecule in self._scaled_molecules[lambda_i].keys():
            for section in self._scaled_molecules[lambda_i][molecule].keys():
               if section=='header':
                  lines+=self._scaled_molecules[lambda_i][molecule][section]+['\n']
               else:
                  lines+=[f"[ {section} ]\n"]+self._scaled_molecules[lambda_i][molecule][section]+['\n']
         lines+=self._sections['system']+['\n']
         lines+=self._sections['molecules']+['\n']
         self._sections_out[lambda_i]=lines
   
   def _get_scale_nonbonded(self, lambda_on=True):  
      if lambda_on: 
         atomtypes_new = {i:[] for i in self.lambdai}
         nonbonded_new = {i:[] for i in self.lambdai}
         
         for _at in self._atomtypes:
            for lambdai in self.lambdai:
               e_scaled = float32(lambdai) * float32(_at[-1])
               stringout = f'{_at[0]:<12} {_at[1]:<6} {_at[2]:>6} {_at[3]:>8}{_at[4]:^5}{_at[5]:>11}{_at[6]:>13}\n'
               atomtypes_new[lambdai].append(stringout)
               stringout = f'{"s"+_at[0]:<12} {_at[0]:<6} {_at[2]:>6} {_at[3]:>8}{_at[4]:^5}{_at[5]:>11}{e_scaled:>13.6e}\n'
               atomtypes_new[lambdai].append(stringout)
               
         for nbp in self._sections['nonbond_params']:
            nbp_ = nbp.split()
            if '[' in nbp_[0] or '\n' in nbp[:3]:
               continue
            elif ';' in nbp_[0]:
               continue
            elif len(nbp.split()) == 5:
               for lambdai in self.lambdai:
                  nonbonded_new[lambdai].append(nbp)
                  eps = float32(nbp_[4])
                  eps_scaled = float32(lambdai) * eps
                  stringout = f'{"s"+nbp_[0]:>5} {"s"+nbp_[1]:>4} {nbp_[2]:>5} {nbp_[3]:>10} {eps_scaled:>8.6e}\n'
                  nonbonded_new[lambdai].append(stringout)
                  eps_scaled = (float32(lambdai)**0.5) * eps
                  stringout = f'{nbp_[0]:>5} {"s"+nbp_[1]:>4} {nbp_[2]:>5} {nbp_[3]:>10} {eps_scaled:>8.6e}\n'
                  nonbonded_new[lambdai].append(stringout)
                  stringout = f'{"s"+nbp_[0]:>5} {nbp_[1]:>4} {nbp_[2]:>5} {nbp_[3]:>10} {eps_scaled:>8.6e}\n'
                  nonbonded_new[lambdai].append(stringout)
            else: print(f'here is the problem line: {nbp} {nbp_[0]}')
         self._scaled_nonbonded = nonbonded_new
         self._scaled_atomtypes = atomtypes_new
      else:
         self._scaled_nonbonded = {i:deepcopy(self._sections['nonbond_params'][1:]) for i in self.lambdai}
         self._scaled_atomtypes = {i:deepcopy(self._sections['atomtypes'][1:]) for i in self.lambdai}
      pass
   
   def _kappa_atoms(self, **kwargs):
      kappa_molecules = kwargs.get('kappa_molecules', False)
      kappa_atom_names = kwargs.get('kappa_atom_names', False)
      if kappa_atom_names and isinstance(kappa_atom_names, list) and \
         all([ True if isinstance(i, str) else False for i in kappa_atom_names]):
            return array(kappa_atom_names)
      if not (kappa_molecules and isinstance(kappa_molecules, list) and \
         all([ True if isinstance(i, int) else False for i in kappa_molecules])):
         self.show_molecule_names()
         atom_names = []
         queryK = ['How many molecules are getting scaled with kappa (commonly just the solvent): ', "int"]
         queryKsol = ['Add molecule index:', "int"]
         kappa_molecules = molecule_selections(queryK, queryKsol)
      for kappa_molecule_i in kappa_molecules:
         queryKmol = [f'How many atoms of molecule {kappa_molecule_i} ({self.molecules[kappa_molecule_i]}) will be selected:', "int"]
         queryKmol_atoms = [f'Add atom index:', "int"]
         self.show_molecule_atomtypes(kappa_molecule_i)
         for atom_idx in molecule_selections(queryKmol, queryKmol_atoms): 
            atom_names.append(get_molecule_atomtypes(kappa_molecule_i)[atom_idx])
      
      return array(atom_names)
         
   def _cold_atoms(self):
      cold = []
      for i in self.molecules.keys():
         if i not in self.hot_molecules:
            for atoms in get_molecule_atomtypes(self._sections, i):
               cold.append(atoms)
      return cold
   
   def _gather_atomtypes_scaled_idx(self, atomtypes):
      idx = {}
      for num, line in enumerate(self._scaled_atomtypes.__getitem__(1)):
         for atomtype in atomtypes:
            if atomtype in line and 's'+atomtype not in line:
               idx[atomtype] = num
      return idx
   
   def _gather_atomtypes_unscaled_lines(self, atomtypes):
      atomtypes_lines = {}
      for line in self._atomtypes:
         if any(item in line for item in atomtypes):
            atomtypes_lines[line[0]] = line
      return atomtypes_lines
   
   def _replace_scaled_atomtype(self, num, line, lambdai):
      eps = float32(line[-1]) * float32(self.kappa[lambdai])**2
      stringout = f'{line[0]:<12} {line[1]:<6} {line[2]:>6} {line[3]:>8}{line[4]:^5}{line[5]:>11}{eps:>13.6e}\n'
      self._scaled_atomtypes[lambdai][num] = stringout
      pass
   
   def _append_nonbonded_kappa_fix(self, kappa_line, cold_lines):
      append_nonbonded_lines = []
      for cold_line in cold_lines: 
         line = generate_nonbonded(kappa_line,cold_line)
         append_nonbonded_lines.append(line)
      for lambda_i in self.lambdai: self._scaled_nonbonded[lambda_i] += append_nonbonded_lines
      pass
   
   def _generate_nonbonded_kappa_fix(self, **kwargs):
      cold_atoms = self._cold_atoms()
      kappa_atoms = self._kappa_atoms(**kwargs) if not self.kappa_atoms else self.kappa_atoms
      assert isinstance(kappa_atoms, (list,ndarray))
      kappa_atomtypes_lines = self._gather_atomtypes_unscaled_lines(kappa_atoms)
      cold_atomtypes_lines = self._gather_atomtypes_unscaled_lines(cold_atoms)
      kappa_atomtypes_scaled_idx = self._gather_atomtypes_scaled_idx(kappa_atoms)
      for kappa_atom in kappa_atoms:
         self._append_nonbonded_kappa_fix( kappa_atomtypes_lines[kappa_atom], cold_atomtypes_lines.values())
         for lambda_i in self.lambdai: 
            self._replace_scaled_atomtype(kappa_atomtypes_scaled_idx[kappa_atom], kappa_atomtypes_lines[kappa_atom], lambda_i)
      pass
   
   def show_molecule_atomtypes(self, molecule:int=0):
      for i in enumerate(get_molecule_atomtypes(self._sections, molecule)): print('{}: {}'.format(*i))
      pass
   
   def _ordered_molecules(self):
      try:
         molecules = [" ".join(i.split()[:-1]) for i in self._sections['molecules'] if '[' not in i and ';' not in i.split()[0]] 
         self.molecules = {num: i for num,i in enumerate(molecules)}
         self.nmol = max(self.molecules.keys())
      except:
         print("Do you have molecules?")
      pass
   
   def show_molecule_names(self):
      try:
         print("\n".join([f'{key}: {mol}' for key,mol in zip(self.molecules.keys(),self.molecules.values())]))
      except:
         print("Do you have molecules?")
      pass
   
   def _moleculetype_sub(self,linestart:int):
      first_round = True
      output = []
      for line in self.readfile[linestart:]:
         if "[" not in line and ';' not in line[:3]:
            output.append(line)
         elif ";" == line[0]:
            continue
         elif "[" in line and first_round!=True: 
            break
         first_round = False
      return output 
   
   def _identify_moltype_sections(self,trunks:int):
      section_start = []
      for i, line in enumerate(self.readfile[trunks:]):
         if '[' in line and 'moleculetype' not in line and 'system' not in line:
            section_start.append(i+trunks)
         elif i != 0 and 'moleculetype' in line or 'system' in line:
            break
      return section_start

   def _parse_moleculetypes(self,trunks:int):
      output = {}
      sections = self._identify_moltype_sections(trunks)
      output['header'] = self.readfile[trunks:trunks+3]
      for section in sections:
         output[self.readfile[section].split()[1]] = []
      for section in sections:
         section_ = self.readfile[section].split()[1]
         output[section_] += self._moleculetype_sub(section)
      return output

   def _gather_param_sections(self):
      for section in self.hard_order_sections[:-1]:
         is_select = [ i for i, line in enumerate(self.readfile) if section in line ]
         in_select = [f' [ {section} ] \n']
         for i in is_select:
            in_select += _parse_section(self.readfile,i)
         self._sections[section]=in_select
      section = self.hard_order_sections[-1]
      is_select = [ i for i, line in enumerate(self.readfile) if section in line ] 
      moltype_dict = {} 
      for i in range(len(is_select)):
         moltype_dict[i] = self._parse_moleculetypes(is_select[i])
      self._sections[section] = moltype_dict
      pass
   
   def _request_hot_molecules(self,**kwargs):
      self.hot_molecules = kwargs.get('hot_molecules', None)
      kappa_method = True if 'kappa' in self._methods[kwargs.get('method', None)] else False
      
      if not self.hot_molecules and not isinstance(self.hot_molecules, list):
         queryN = ["How many molecules are you going to temper?", "int"]
         queryM = ["Add molecule index:", "int"]
         print("At least 1 molecule must be selected as the hot molecule.")
         if kappa_method:
            print("For methods with kappa, select the solute(s) to scale nonbonded water interactions.")
            print("All non-selected molecules will have interactions reset to processed.top ff terms.")
         self.hot_molecules = molecule_selections(queryN, queryM)
      pass
   
   def _request_kappa_params(self, **kwargs):
      self.kappa_low_temp = kwargs.get('kappa_low_temp', self.kappa_low_temp)
      self.kappa = kwargs.get('kappa_max', None)
      if not isinstance(self.kappa, dict):
         self.kappa = float(input('Enter the inputs (default: 1.06): ').strip() or "1.06")
      pass
   
   def _method_rest2(self, **kwargs):
      if kwargs.get('verbose',False): print(f'Running rest2 scaling method')
      temps = kwargs.get('temps', None)
      temps_opt = kwargs.get('temps_opt',None)
      if temps_opt is not None:
         if len(temps_opt) == self.nreps:
            self.lambdai = self.templadder
         else: raise Exception(f'Optional temp ladder {len(temps)} does not match the number of replicas {self.nreps}')
      elif self._processed:
         if temps:
            self.templadder = generate_temperature_ladder(self.nreps, *temps)
            self.lambdai = self.templadder
         else:
            raise TypeError('must provide kwargs entry \"temps\" (2,) when running calculations in sucession.')
      elif temps:
         self.templadder = generate_temperature_ladder(self.nreps,*temps)
         self.lambdai = self.templadder
      elif self.templadder:
         assert len(self.templadder) == self.nreps, f'Number of replicas {self.nreps} and temperature ladder {len(self.templadder)}'
         self.lambdai = self.templadder
      assert len(self.templadder) == self.nreps, f'Number of replicas {self.nreps} and temperature ladder {len(self.templadder)}'
      self._request_hot_molecules(**kwargs)
      self._get_scale_nonbonded()
      self._get_scale_dehedrals()
      self._get_scaled_molecules()
      self._populate_out()
      self._processed = True
      pass
   
   def _method_solvent_scaled(self,**kwargs):
      if kwargs.get('verbose',False): print(f'Running solvent scaling method')
      self.nreps = kwargs.get('nreps', self.nreps)
      self._lambdai = arange(1,self.nreps+1)
      self._request_hot_molecules(**kwargs)
      self._get_scale_nonbonded(lambda_on=False)
      self.kappa_low_temp = kwargs.get('kappa_low_temp', self.kappa_low_temp)
      self._request_kappa_params(**kwargs)
      self._generate_nonbonded_kappa_fix(**kwargs)
      self._get_scale_dehedrals(lambda_on=False)
      self._get_scaled_molecules(lambda_on=False)
      self._populate_out()
      self._processed = True
      pass
   
   def _method_ssrest3(self, **kwargs):
      if kwargs.get('verbose',False): print(f'Running ssrest3 scaling method')
      temps = kwargs.get('temps', None) if not any(self.templadder) else None
      if temps: self.templadder = generate_temperature_ladder(self.nreps,*temps)
      temps_opt = kwargs.get('temps_opt',None)
      if temps_opt is not None:
         if len(temps_opt) == self.nreps:
            self.lambdai = self.templadder
         else: raise Exception(f'Optional temp ladder {len(temps)} does not match the number of replicas {self.nreps}')
      elif temps: 
         self.templadder = generate_temperature_ladder(self.nreps,*temps)
         self.lambdai = self.templadder
      if self._processed:
         self.lambdai = self.templadder
      self._request_hot_molecules(**kwargs)
      self._get_scale_nonbonded()

      self._request_kappa_params(**kwargs)
      self._generate_nonbonded_kappa_fix(**kwargs)
      self._get_scale_dehedrals()
      self._get_scaled_molecules()
      self._populate_out()
      self._processed = True
      pass
  
   def run(self, **kwargs):
      self.nreps = kwargs.get('nreps', self.nreps)
      scaling_method = kwargs.get('method', None)
      if kwargs.get('verbose',False):
         print(scaling_method)
      if scaling_method == None:
         raise Exception(f'No scaling method selected.\n Options used: {kwargs}')
      elif scaling_method == 'solvent_scaled':
         self._method_solvent_scaled(**kwargs)
      elif scaling_method == 'rest2':
         self._method_rest2(**kwargs)
      elif scaling_method == 'ssrest3':
         self._method_ssrest3(**kwargs)
      else:
         raise Exception(f'\"{scaling_method}\" is not an implemented method.')

      outfile = kwargs.get('outfile', None)
      if outfile:
         for lambda_i in self.lambdai:
            if 'kappa' in self._methods[scaling_method]:
               prop = {'outfile': outfile,
                        'filepath': kwargs.get('filepath','./'),
                        'lambda_i': lambda_i,
                        'kappa_i': self.kappa[lambda_i],
                        'lines': self._sections_out[lambda_i]}
            else:
               prop = {'outfile': outfile,
                        'filepath': kwargs.get('filepath','./'),
                        'lambda_i': lambda_i,
                        'lines': self._sections_out[lambda_i]}
            topology_writer(**prop)
      pass 
   
if __name__=="__main__":
   import argparse
   parser = argparse.ArgumentParser(
                    prog='repex_topology_parser',
                    description='Parse GROMACS processed topology file and perform desired scaling.',
                    epilog='Happy scaling!')
   parser.add_argument('-p','--topol', default='processed.top', type=str, help='Processed Topology (default: %(default)s)')
   parser.add_argument('--show-molecules', action='store_true', help='Show molecule indices found in processed.top')
   parser.add_argument('--show-molecule-atomtypes', default=None, type=int, help='Show molecule atomtypes found in processed.top')
   parser.add_argument('--verbose', action='store_true')
   parser.add_argument('-O','--outfile', default='topol', type=str, help='Scaled topology base name (default: %(default)s)')
   parser.add_argument('-P','--output', default='./', type=str, help='Scaled topology output PATH (default: %(default)s)')
   parser.add_argument('-m','--method', default='rest2', type=str, help='Scaled topology method [ rest2 | ssrest3 ] (default: %(default)s)')
   parser.add_argument('-H','--hot-molecules', default=[0], nargs='+', type=int, help='Space delimited entry of hot molecule selection (default: %(default)s)')
   parser.add_argument('-n','--nreps', default=10, type=int, help='Number of replicas (default: %(default)s)')
   parser.add_argument('-k', '--kappa-max', default=1.06, type=float, help='Max Kappa Scaling (default: %(default)s)')  
   parser.add_argument('--kappa-low-temp', default=300, type=float, help='Replicas at or above this temperature will have active solvent scaling (default:%(default)s) ')  
   # parser.add_argument('-t','--temps-opt', default=[], nargs='+', type=float, help='Space delimited entry of user selected temperature ladder (default: %(default)s)')
   parser.add_argument('--kappa-atomtypes', default=['OW'], nargs='+', type=str, help='Atomtypes to apply solvent scaling (default:%(default)s) ') 
   parser.add_argument('-T', '--tmin', default=300, type=float, help='Base Temperature (default: %(default)s)')
   parser.add_argument('-M', '--tmax', default=500, type=float, help='Max Temperature (default: %(default)s)')
   
   args = parser.parse_args()
   assert args.tmin < args.tmax, 'Base replica temperature must be lower than hot replica!'
   options = {'ifile':args.topol, 
              'outfile':args.outfile, 
              'filepath':args.output, 
              'method':args.method, 
              'hot_molecules':args.hot_molecules, 
              'nreps':args.nreps, 
              'kappa_max':args.kappa_max, 
              'temps':[args.tmin,args.tmax] , 
              'kappa_low_temp': args.kappa_low_temp, 
              'kappa_atom_names':args.kappa_atomtypes,
              'verbose':args.verbose,
            #   'temps_opt':args.temps_opt
              }
   
   conversion_class = topo2rest(**options)
   if args.show_molecules:
      conversion_class.show_molecule_names()
   elif args.show_molecule_atomtypes:
      conversion_class.show_molecule_atomtypes(args.show_molecule_atomtypes)
   else:
      conversion_class.run(**options)
   