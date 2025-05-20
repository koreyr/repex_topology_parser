#!/bin/env python3



from numpy import exp, log, arange, array, float32, unique, sqrt, where, ndarray, isin
import pandas as pd
from io import StringIO
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

def query_typed(query, qtype):
   return eval(qtype)(input(query))

def molecule_selections(queryN:list, queryM:list):
   response = query_typed(*queryN)
   molecule_sel = []
   for i in range(response):
      molecule_sel.append(query_typed(*queryM))
   return molecule_sel

def topology_writer(**kwargs):
   ofile = kwargs.get('outfile', 'topol')
   filepath = kwargs.get('filepath', './')
   lambda_i = kwargs.get('lambda_i', None)
   kappa_i = kwargs.get('kappa_i', None)
   lines = kwargs.get('lines', None)
   if any(lines):
      fullpath = f'{filepath}{ofile}'
      for additional in [lambda_i, kappa_i]: fullpath = f'{fullpath}-{additional:0.3f}' if additional else fullpath
      fullpath += ".top"
      with open(fullpath,'w') as file:
         file.writelines(lines)
      print(f'Saved {fullpath}')
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
      self._kappa = None
      self.kappa_atoms = None
      self.hot_molecules = None
      self.nmol = None
      self.molecules = None
      self.templadder = generate_temperature_ladder(nreps, *temps)
      self.lambdai = compute_lambda(self.templadder)
      self._sections = {}
      self._sections_out = {}
      self.scaled_dihedrals = {}
      self.scaled_dihedral_types = {}
      self._scaled_atomtypes = {}
      self._scaled_charges={}
      self.scaled_nonbonded={}
      self._molecule_atoms = {}
      self.kappa_low_temp = 330
      self.kappa = 1.0
      
      self.hard_order_sections = [ "defaults", "atomtypes", "nonbond_params", "bondtypes", \
                                   "constrainttypes", "angletypes", "dihedraltypes","system", "molecules", "moleculetype"]
      with open(ifile) as topo:
         self.readfile = topo.readlines()
      self._gather_param_sections()
      self._atomtypes = [i[:i.find(';')].split() \
                               for i in self._sections['atomtypes'] if ';' not in i[:3] if '[' not in i[:3] \
                               if '\n' not in i[:3]]
      self.molecule_atoms = True
      self.show_molecule_names()
       
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
         print('made it')
         atdict = {}
         for nmol in range(len(self._sections['moleculetype'])):
            atdict[nmol] = {int(i.split()[0]): i.split()[1] for i in self._sections['moleculetype'][nmol]['atoms'] if len(i.split()) != 0 and ';' not in i[:5]}
         self._molecule_atoms = atdict
      pass
      
   def _get_scaled_charges(self):
      for l in self.lambdai:
         self._scaled_charges[l]=deepcopy(self._sections['moleculetype'])
         for hot in self.hot_molecules:
            a=[i.split() if len(i.split()) == 8 else i.split()[:-3] for i in self._sections['moleculetype'][hot]['atoms'] if len(i.split()) != 0]
            c=[]
            for i in a:
               s_charge=float(i[6])*sqrt(l) # scalling with sqrt(lambda)
               stringout = f'{i[0]:>6} {"s"+i[1]:>10} {i[2]:>6} {i[3]:>6} {i[4]:>6} {i[5]:>6} {s_charge:>10.6e} {float(i[7]):>10.3f}\n'
               c.append(stringout)   
            self._scaled_charges[l][hot]['atoms'] = c
            self._scaled_charges[l][hot]['dihedrals'] = deepcopy(self.scaled_dihedrals[hot][l])

   def _populate_out(self):
      for lambda_i in self.lambdai:
         lines=self._sections['defaults']+['\n']
         lines+=['[ atomtypes ]\n']+self._scaled_atomtypes[lambda_i]+['\n']
         lines+=['[ nonbond_params ]\n']+self.scaled_nonbonded[lambda_i]+['\n']
         lines+=self._sections['bondtypes']+['\n']
         lines+=self._sections['constrainttypes']+['\n']
         lines+=self._sections['angletypes']+['\n']
         lines+=self._sections['dihedraltypes']+['\n']
         for i in self._scaled_charges[lambda_i].keys():
            for j in self._scaled_charges[lambda_i][i].keys():
               if j=='header':
                  lines+=self._scaled_charges[lambda_i][i][j]+['\n']
               else:
                  lines+=[f"[ {j} ]\n"]+self._scaled_charges[lambda_i][i][j]+['\n']
         lines+=self._sections['system']+['\n']
         lines+=self._sections['molecules']+['\n']
         self._sections_out[lambda_i]=lines
   
   def _parse_section(self,trunks):
      first_round = True
      output = []
      for line in self.readfile[trunks:]:
         if "[" not in line and first_round!=True and len(line.split()) != 0:
            output.append(line)
         elif ";" == line[0]: 
            continue
         elif "[" in line and first_round!=True: 
            break
         first_round = False
      return output

   def _get_molecule_atomtypes(self,molecule:int=0):
      b=[]
      for i in self._sections['moleculetype'][molecule]['atoms']:
         if len(i.split())>1 and ';' not in i.split()[0] and '[' not in i :
            b.append(i.split()[:2])
      dataset = array(b,dtype=object)
      atomtypes = dataset[:,1]
      return unique(atomtypes)
   
   def _get_scale_nonbonded(self):   
      atomtypes_new = {}
      nonbonded_new = {}
      for lambdai in self.lambdai:
         at_out = []
         for _at in self._atomtypes:
            e_scaled = float32(lambdai) * float32(_at[-1])
            stringout = f'{_at[0]:<12} {_at[1]:<6} {_at[2]:>6} {_at[3]:>8}{_at[4]:^5}{_at[5]:>11}{_at[6]:>13}\n'
            at_out.append(stringout)
            stringout = f'{"s"+_at[0]:<12} {_at[0]:<6} {_at[2]:>6} {_at[3]:>8}{_at[4]:^5}{_at[5]:>11}{e_scaled:>13.6e}\n'
            at_out.append(stringout)
         atomtypes_new[lambdai] = at_out
         nbp_out = []
         for nbp in self._sections['nonbond_params']:
            nbp_ = nbp.split()
            if '[' in nbp_[0] or '\n' in nbp[:3]:
               continue
            elif ';' in nbp_[0]:
               continue
            elif len(nbp.split()) == 5:
               nbp_out.append(nbp)
               eps = float32(nbp_[4])
               eps_scaled = float32(lambdai) * eps
               stringout = f'{"s"+nbp_[0]:>5} {"s"+nbp_[1]:>4} {nbp_[2]:>5} {nbp_[3]:>10} {eps_scaled:>8.6e}\n'
               nbp_out.append(stringout)
               eps_scaled = (float32(lambdai)**0.5) * eps
               stringout = f'{nbp_[0]:>5} {"s"+nbp_[1]:>4} {nbp_[2]:>5} {nbp_[3]:>10} {eps_scaled:>8.6e}\n'
               nbp_out.append(stringout)
               stringout = f'{"s"+nbp_[0]:>5} {nbp_[1]:>4} {nbp_[2]:>5} {nbp_[3]:>10} {eps_scaled:>8.6e}\n'
               nbp_out.append(stringout)
            else: print(f'here is the problem line: {nbp} {nbp_[0]}')
         nonbonded_new[lambdai] = nbp_out
      self.scaled_nonbonded = nonbonded_new
      self._scaled_atomtypes = atomtypes_new
      pass
   
   def _kappa_atoms(self):
      self.show_molecule_names()
      atom_names = []
      queryK = ['How many molecules are getting scaled with kappa (commonly just the solvent): ', "int"]
      queryKsol = ['Add molecule index:', "int"]
      for kappa_molecule_i in molecule_selections(queryK, queryKsol):
         queryKmol = [f'How many atoms of molecule {kappa_molecule_i} ({self.molecules[kappa_molecule_i]}) will be selected:', "int"]
         queryKmol_atoms = [f'Add atom index:', "int"]
         self.show_molecule_atomtypes(kappa_molecule_i)
         for atom_idx in molecule_selections(queryKmol, queryKmol_atoms): 
            atom_names.append(self._get_molecule_atomtypes(kappa_molecule_i)[atom_idx])
      return array(atom_names)
         
   def _cold_atoms(self):
      cold = []
      for i in self.molecules.keys():
         if i not in self.hot_molecules:
            for atoms in self._get_molecule_atomtypes(i):
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
   
   def _replace_scaled_atomtype(self, num, line, kappa):
      eps = float32(line[-1]) * float32(kappa)**2
      stringout = f'{line[0]:<12} {line[1]:<6} {line[2]:>6} {line[3]:>8}{line[4]:^5}{line[5]:>11}{eps:>13.6e}\n'
      self._scaled_atomtypes[num] = stringout
      pass
   
   def _generate_nonbonded(self, at1, at2):
      e1 = float32(at1[-1])
      e2 = float32(at2[-1])
      s1 = float32(at1[-2])
      s2 = float32(at2[-2])
      funct = str(1)
      eij, sij = compute_se(e1,e2,s1,s2)
      a1, a2 = at1[0], at2[0]
      return f'{a1:>5} {a2:>4} {funct:^5} {sij:>10.4f} {eij:>8.4f}\n'
   
   def _append_nonbonded_kappa_fix(self, no_kappa_line, cold_lines):
      append_nonbonded_lines = []
      for cold_line in cold_lines: 
         line = self._generate_nonbonded(no_kappa_line,cold_line)
         append_nonbonded_lines.append(line)
      for lambda_i in self.lambdai: self.scaled_nonbonded[lambda_i] += append_nonbonded_lines
      pass
   
   def _generate_nonbonded_kappa_fix(self):
      cold_atoms = self._cold_atoms()
      kappa_atoms = self._kappa_atoms() if not self.kappa_atoms else self.kappa_atoms
      assert isinstance(kappa_atoms, (list,ndarray))
      kappa_atomtypes_lines = self._gather_atomtypes_unscaled_lines(kappa_atoms)
      cold_atomtypes_lines = self._gather_atomtypes_unscaled_lines(cold_atoms)
      kappa_atomtypes_scaled_idx = self._gather_atomtypes_scaled_idx(kappa_atoms)
      for kappa_atom in kappa_atoms:
         self._append_nonbonded_kappa_fix( kappa_atomtypes_lines[kappa_atom], cold_atomtypes_lines.values())
         for lambda_i in self.lambdai: 
            self._replace_scaled_atomtype(kappa_atomtypes_scaled_idx[kappa_atom], kappa_atomtypes_lines[kappa_atom], self.kappa[lambda_i])
      pass
   
   def show_molecule_atomtypes(self, molecule:int=0):
      for i in enumerate(self._get_molecule_atomtypes(molecule)): print('{}: {}'.format(*i))
      pass
   
   def show_molecule_names(self):
      try:
         molecules = [" ".join(i.split()[:-1]) for i in self._sections['molecules'] if '[' not in i and ';' not in i.split()[0]] 
         molecules = {num: i for num,i in enumerate(molecules)}
         print("\n".join([f'{key}: {mol}' for key,mol in zip(molecules.keys(),molecules.values())]))
         if self.molecules == None: self.molecules = molecules
         if self.nmol == None: self.nmol = max(self.molecules.keys())
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
   
   def _get_scale_dehedrals_(self):

      dihedrals_new = {}
      dihedral_types_new = {}
      dihedral_types = [" ".join(i.split()[:-2]) if i.split()[-2] == ';' else i.split()[:-1] if i.split()[-1] == ';' \
                               else i.split() for i in self._sections['dihedraltypes'] if ';' not in i[:3] if '[' not in i[:3] \
                               if '\n' not in i[:3]]
      
      for lambdai in self.lambdai:
         dih_types_new = []
         for dihedraltype in dihedral_types:
            dtls_ = dihedraltype.split()
            Kscaled = float32(dtls_[6])*lambdai
            
            if len(dtls_) == 8:
               check_atom_X = isin(array(dtls_[:4]), array(['X']))
               if sum(check_atom_X) == 0:
                  stringout = f'{"s"+dtls_[0]:<5} {"s"+dtls_[1]:<5} {"s"+dtls_[2]:<5} {"s"+dtls_[3]:<5} {dtls_[4]:^9}{dtls_[5]:<10}{Kscaled:<10.5f}{dtls_[7]}\n'
                  dih_types_new.append(stringout)
               elif dtls_[0] == 'X':
                  if sum(check_atom_X) == 1:
                     stringout = f'{dtls_[0]:<5} {"s"+dtls_[1]:<5} {"s"+dtls_[2]:<5} {"s"+dtls_[3]:<5} {dtls_[4]:^9}{dtls_[5]:<10}{Kscaled:<10.5f}{dtls_[7]}\n' 
                     dih_types_new.append(stringout)
                  elif dtls_[3] == 'X':
                     stringout = f'{dtls_[0]:<5} {"s"+dtls_[1]:<5} {"s"+dtls_[2]:<5} {dtls_[3]:<5} {dtls_[4]:^9}{dtls_[5]:<10}{Kscaled:<10.5f}{dtls_[7]}\n' 
                     dih_types_new.append(stringout)
                  elif dtls_[1] == 'X':
                     stringout = f'{dtls_[0]:<5} {dtls_[1]:<5} {"s"+dtls_[2]:<5} {"s"+dtls_[3]:<5} {dtls_[4]:^9}{dtls_[5]:<10}{Kscaled:<10.5f}{dtls_[7]}\n' 
                     dih_types_new.append(stringout)
                  else: print(f'In dihedraltype: something slipped through\n{dihedraltype}')
                  
               else: print(f'warning: parameter not found for {dtls_[:4]}')
            else: print('warning: dihedraltype not processed:\n '+dihedraltype)
         dihedral_types_new[lambdai] = dih_types_new
      self.scaled_dihedral_types = dihedral_types_new
      
      for hot in self.hot_molecules:
         dihedrals = self._sections['moleculetype'][hot]['dihedrals']
         atn2t = self.molecule_atoms[hot]
         for lambdai in self.lambdai:
            df_dht = pd.read_csv(StringIO(''.join(self.scaled_dihedral_types[lambdai])),names=['i','j','k','l','func','angle','K','mult'], sep='\s+')
            # attempt to convert code to avoid issue, something is wrong
            # df_dht = pd.DataFrame([ i.split() for i in self.scaled_dihedral_types[lambdai]], columns=['i','j','k','l','func','angle','K','mult'])
            dih_new = []
            warn_once = True
            for dihedral in dihedrals:
               dls_ = dihedral[:dihedral.find(';')].split()
               if len(dls_) == 5:
                  i, j, k, l = [atn2t[int(atomnum)] for atomnum in dls_[:4]]
                  func = int(dls_[4])
                  dht_series = df_dht[(df_dht['i']=='s'+i) * (df_dht['j']=='s'+j) * (df_dht['k']=='s'+k) * (df_dht['l']=='s'+l) * (df_dht['func']==func)]
                  dht_series = pd.concat( [dht_series,df_dht[(df_dht['i']=='s'+l) * (df_dht['j']=='s'+k) * (df_dht['k']=='s'+j) * (df_dht['l']=='s'+i) * (df_dht['func']==func)]],ignore_index=True)
                  if func==4 and dht_series.empty:
                     dht_series = pd.concat( [dht_series,df_dht[(df_dht['i']=='X') * (df_dht['j']=='X') * (df_dht['k']=='s'+k) * (df_dht['l']=='s'+l) * (df_dht['func']==func)]],ignore_index=True)
                     dht_series = pd.concat( [dht_series,df_dht[(df_dht['i']=='X') * (df_dht['j']=='X') * (df_dht['k']=='s'+j) * (df_dht['l']=='s'+i) * (df_dht['func']==func)]],ignore_index=True)
                     dht_series = pd.concat( [dht_series,df_dht[(df_dht['i']=='X') * (df_dht['j']=='s'+j) * (df_dht['k']=='s'+k) * (df_dht['l']=='s'+l) * (df_dht['func']==func)]],ignore_index=True)
                     dht_series = pd.concat( [dht_series,df_dht[(df_dht['i']=='X') * (df_dht['j']=='s'+k) * (df_dht['k']=='s'+j) * (df_dht['l']=='s'+i) * (df_dht['func']==func)]],ignore_index=True)
                  if func==9 and dht_series.empty:
                     if j==k:
                        dht_series = pd.concat( [dht_series,df_dht[(df_dht['i']=='X') * (df_dht['j']=='s'+k) * (df_dht['k']=='s'+j) * (df_dht['l']=='X') * (df_dht['func']==func)]],ignore_index=True)
                     else:
                        dht_series = pd.concat( [dht_series,df_dht[(df_dht['i']=='X') * (df_dht['j']=='s'+j) * (df_dht['k']=='s'+k) * (df_dht['l']=='X') * (df_dht['func']==func)]],ignore_index=True)
                        dht_series = pd.concat( [dht_series,df_dht[(df_dht['i']=='X') * (df_dht['j']=='s'+k) * (df_dht['k']=='s'+j) * (df_dht['l']=='X') * (df_dht['func']==func)]],ignore_index=True)
                  dht_series = dht_series.drop_duplicates()
                  dht_series.reset_index()
                  if dht_series.empty and warn_once:
                     warn_once = False
                     print("Warning: Dihedral not found!")
                     print(i,j,k,l,func,*dls_[:4])
                     print(dht_series)
                  for _, dht_ in dht_series.iterrows():
                     stringout = f'{dls_[0]:<5} {dls_[1]:>5} {dls_[2]:>5} {dls_[3]:>5} {dht_["func"]:^9}{dht_["angle"]:<10}{dht_["K"]:<10.5f}{dht_["mult"]}\n'
                     dih_new.append(stringout)
               elif len(dls_) == 8:
                  Kscaled = float32(dls_[6])*lambdai
                  stringout = f'{dls_[0]:<5} {dls_[1]:>5} {dls_[2]:>5} {dls_[3]:>5} {dls_[4]:^9}{dls_[5]:<10}{Kscaled:<10.5f}{dls_[7]}\n'
                  dih_new.append(stringout)
               elif len(dls_) == 0:
                  continue
               else: 
                  print(f'Warning: incorrect parsing of dihedrals section\n expecting 5 or 8 columns found {len(dls_)}\n'+dihedral)
            dihedrals_new[lambdai] = deepcopy(dih_new)
         self.scaled_dihedrals[hot] = deepcopy(dihedrals_new) 
      pass

   def _gather_param_sections(self):
      for section in self.hard_order_sections[:-1]:
         is_select = [ i for i, line in enumerate(self.readfile) if section in line ]
         in_select = [f' [ {section} ] \n']
         for i in is_select:
            in_select += self._parse_section(i)
         self._sections[section]=in_select
      section = self.hard_order_sections[-1]
      is_select = [ i for i, line in enumerate(self.readfile) if section in line ] 
      moltype_dict = {} 
      for i in range(len(is_select)):
         moltype_dict[i] = self._parse_moleculetypes(is_select[i])
      self._sections[section] = moltype_dict
      pass
  
   def run(self, **kwargs):
      self.nreps = kwargs.get('nreps', self.nreps)
      temps = kwargs.get('temps', None) if not any(self.templadder) else None
      if temps: self.templadder = generate_temperature_ladder(temps)
      self.hot_molecules = kwargs.get('hot_molecules', None)
      if not self.hot_molecules and not isinstance(self.hot_molecules, list):
         queryN = ["How many molecules are you going to temper?", "int"]
         queryM = ["Add molecule index:", "int"]
         hot_molecules = kwargs.get('hot_molecules', None )
         self.hot_molecules = hot_molecules if not hot_molecules and isinstance(hot_molecules, list) \
                                            else molecule_selections(queryN, queryM)
      self._get_scale_nonbonded()
      if kwargs.get('apply_kappa', False):
         self.kappa_low_temp = kwargs.get('kappa_low_temp', self.kappa_low_temp)
         self.kappa = kwargs.get('kappa_max', None)
         if not isinstance(self.kappa, dict):
             self.kappa = float(input('Enter the inputs (default: 1.06): ').strip() or "1.06")
         self._generate_nonbonded_kappa_fix()
      self._get_scale_dehedrals_()
      self._get_scaled_charges()
      self._populate_out()
      outfile = kwargs.get('outfile', None)
      if outfile:
         for lambda_i in self.lambdai:
            topology_writer(\
                        **{'outfile': outfile,
                        'filepath': kwargs.get('filepath','./'),
                        'lambda_i': lambda_i,
                        'kappa_i': self.kappa[lambda_i],
                        'lines': self._sections_out[lambda_i]})
      pass 