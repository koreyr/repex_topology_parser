# Python library to parse and scale processed GROMACS topology files
#### Current version 0.2.2 (support for amber potentials, CHARMM and OPLS tbd)
### Correctly Scales topologies for:
- Replica exchange with solute scaling (REST2)
- Solvent Scaled Replica Exchange with solute scaling (ssREST3)

#### Example files are included for REST2 and ssREST3 test runs
##### Code supports solute scaling of multiple solutes via interactive selection
##### Code supports solvent scaling of multiple solutes and multiple atoms via interactive selection


##### ============================================================
### repex_topology_parser
#### Currently only supports Amber Potential 
#### without CMAP corrections
```math
\begin{align*}
E_{\rm total} =  \sum_{\rm bonds} K_r \big(r - r_{\rm eq}\big)^2& + \sum_{\rm angles} K_\theta \big(\theta - \theta_{\rm eq}\big)^2 + \sum_{\rm dihedrals} \frac{V_n}{2} \Big( 1 + \cos(n\phi - \gamma) \Big) \\
&+ \sum_{i\lt j} \epsilon_{ij} \left( \left(\frac{\sigma_{ij}}{R_{ij}} \right)^{12} - \left( \frac{\sigma_{ij}}{R_{ij}} \right)^6 \right) + \sum_{i\lt j} \frac{q_i q_j}{4 \pi \epsilon_0 R_{ij}}
\end{align*}
```
### TODO
#### - Add CMAP lambda scaling
#### - Extend to CHARMM, OPLS-AA
### Examples for REST2, and ssREST3 (solvent-scaled REST2)
##### ============================================================

##### ============================================================
### Energy Components
$$ E_{\rm total} = E_{\rm protein-protein} + E_{\rm protein-water} + E_{\rm water-water} $$
#### REST2 scaling with $\lambda$
$$ \lambda = \frac{\beta_n}{\beta_0} = \frac{T_0}{T_n} : \beta_n = \frac{1}{k_B T_n} $$
$$ E_{\rm total}^{\lambda_n} = \lambda_n E_{\rm protein-protein} + \sqrt{\lambda_n} E_{\rm protein-water} + E_{\rm water-water} $$
### solvent-scaling 
$$ \kappa_i =  e^{\frac{i}{N - 1} log(\kappa_{max})} ; i \in [0,N-1]$$
$$ E_{\rm total}^{\kappa_n} = E_{\rm protein-protein} + \kappa_n E_{\rm protein-water} + E_{\rm water-water} $$
### ssREST2 (solvent-scaled REST2)
$$ E_{\rm total}^{\lambda_n,\kappa_n} = \lambda_nE_{\rm protein-protein} + \boldsymbol{\kappa_n} \sqrt{\lambda_n} E_{\rm protein-water} + E_{\rm water-water} $$

##### ============================================================


```python
# Load our library
import src.repex_topology_parser as rtp
from pathlib import Path
Path("./example_conversion").mkdir(parents=True, exist_ok=True)
```


```python
# initialize our class providing an input processed.top file
example_scaling = rtp.topo2rest('/home/koreyr/github/repex_topology_parser/tests/test_topologies/apo_sys2.top')
```


```python
### Let us display our molecules contained in our topology
example_scaling.molecules
```




    {0: 'Protein_chain_A', 1: 'SOL', 2: 'NA', 3: 'CL'}



### Here we perform with one command rest2 scaling. The hot molecule will have dihedrals/charges/LJ 
### parameters scaled by $\sqrt\lambda$. Thus intra-hot molecule interactions will be scaled
### by $\lambda$, while hot molecule - other molecules will be scaled by $\sqrt\lambda$


```python
# Our first example is performing solute scaling (REST2) on just the protein
REST2 = { 'hot_molecules':[0], # Select the molecule(s) you desire to scale, as a list
            'nreps':20, # define the number of replicas you desire
            'outfile':'topol_rest2_apo_sys2', # provide a prefix name for your topology
                                     # default = 'topol'
            'filepath':'./tests/test_python/', # define the directory you wish to write your scaled topologies
                                      # default='./'
            'method':'rest2', # method is rest2
            'temps':[300,500], # temperature range 300 to 500, utilized to compute the geometric
                               # temperature ladder
                               # default = [300,500]
            'verbose':True     # be verbose, important if you are performing scaling interactively
                               # e.g. not providing one or more of these options
            }
example_scaling.run(**REST2)
```

    Running rest2 scaling method


### Here we perform ssrest3 scaling with one command, ssrest3 applies additional scaling on the OW atom of our water model, in this case \'OW_tip4pd\' atomtype. 
### First, the input of a hot molecule has the same effect as rest2 where the protein dihedrals/charges/LJ parameters are scaled by $\lambda$.
### Second the LJ $\epsilon$ of water is scaled by $\kappa^2$ tuning the solvation of the hot molecule(s).
### And lastly, all non-hot molecule-water LJ parameters are reset thus avoiding the effects of solvent scaling. 


```python
# Our second example is performing solvent-scaling REST3 (ssREST3) on just the protein
# From the displayed atomtypes contained in our solvent molecule we opt to scale
# 'OW_tip4pd' ('HW' and 'MW' have epsilon = 0.0 so we exclude these from the list)
ssREST3 = { 'hot_molecules':[0], # Select the molecule(s) you desire to scale, as a list
            'nreps':20, # define the number of replicas you desire
            'outfile':'topol_ssrest3_apo_sys2', # provide a prefix name for your topology
                                     # default = 'topol'
            'filepath':'./example_conversion/', # define the directory you wish to write your scaled topologies
                                      # default='./'
            'method':'ssrest3', # method is ssrest3
            'temps':[300,500], # temperature range 300 to 500, utilized to compute the geometric
                               # temperature ladder
                               # default = [300,500]
            'kappa_low_temp' : 300, # At which temperature to activate solvent scaling, maybe useful
                                    # to increse to 330 when using a lower number of replicas so the base replica
                                    # experiences solvation more accurately. For ssREST3 simulations with 
                                    # 16 or more replicas, it is unlikedly changing this value to 330 will have 
                                    # any benefit. 
            'kappa_max' : 1.1, # Maximum kappa value
            'kappa_atom_names' : ['OW_tip4pd'], # List of atomtypes to apply kappa scaling
            'verbose':True     # be verbose, important if you are performing scaling interactively
                               # e.g. not providing one or more of these options
            }
example_scaling.run(**ssREST3)
```

    Running ssrest3 scaling method


##### The most prudent approach is to perform time continuous simulations of the hot or hightemp replica applying different solvent scaling. 
##### In this example the hot replica's effective temperature is 500 K.
##### The conformational sampling of each hightemp replica's continuous simulation should be compared to a continuous unscaled simulation performed with a bath temperature matching the selected effective temperature of the hightemp replica; 200 to 500 ns should do, however longer simulations may be needed.
##### For fully disorderd IDPs: when tested on the 20-residue fragment of $\alpha$-synuclien 200 ns was more than reasonable, a 40-residue fragment of Fused in Sarcoma ~500 ns, and in contrast ~1 $\mu$s for a 70-residue fragment of p27.
##### For this paper, simple metrics were introduced to assess the appropriateness of $\kappa_{max}$ such as radius of gyration and secondary structure population, you are invited to test other metrics that are important to your study. 
#####
##### To minimize the real world time required to sample $\kappa_{max}$ values we opted to select 9 values of $\kappa_{max}$ to test with a $\lambda_{max}$ value of 0.6 corresponding to an effective temperature of 500 K. 
##### On a 10 GPU system this results in:
##### -- One unscaled simulation with a bath temperature coupled to 500 K
##### -- 9 scaled simulations with $\kappa_{max}$ ranging from 1.02 to 1.1
##### Gromacs 2022.5 was compiled with mpi and cuda support. 
##### Simulations were all submitted together using the mdrun flag -multidir:
```bash
gmx_gpu_mpi mdrun -v -deffnm run -mutlidir 500K kappa{1..9}
```

#### The code below generatures the unscaled topology and the 9 values of $\kappa_{max}$:
##### 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09 and 1.1.
#### It is important to check the base replica topology against the processed topology.


```python
from numpy import arange
kappa_scan = rtp.topo2rest('tests/topology_files/test_topo/lig_sys2.top', nreps=2)
for k in arange(1.02,1.1,0.01):
    ssREST3 = { 'hot_molecules':[0], 
                'nreps':2, 
                'outfile':'topol_ssrest3',          
                'filepath':'./test_run/', 
                'method':'ssrest3',
                'temps':[300,500], 
                'kappa_low_temp' : 300, 
                'kappa_max' : k, 
                'kappa_atom_names' : ['OW_tip4pd'], 
                'verbose':True     
                }
    kappa_scan.run(**ssREST3)
```

    Running ssrest3 scaling method
    Running ssrest3 scaling method
    Running ssrest3 scaling method
    Running ssrest3 scaling method
    Running ssrest3 scaling method
    Running ssrest3 scaling method
    Running ssrest3 scaling method
    Running ssrest3 scaling method
    Running ssrest3 scaling method


### Scaling can also be performed from the command line, either you supply the processed topology, the method, the temperature range and if desired $\kappa_{max}$ and enter into an interactive session to select the hot molecule(s) and the solvent to be scaled. 
# Or
### All options can be provided from the command line as shown below. 

```bash
username@computer] python3 src/repex_topology_parser.py --help
usage: repex_topology_parser [-h] [-p TOPOL] [--show-molecules] [--show-molecule-atomtypes SHOW_MOLECULE_ATOMTYPES] [--verbose] [-O OUTFILE] [-P OUTPUT]
                             [-m METHOD] [-H HOT_MOLECULES [HOT_MOLECULES ...]] [-n NREPS] [-k KAPPA_MAX] [--kappa-low-temp KAPPA_LOW_TEMP]
                             [--kappa-atomtypes KAPPA_ATOMTYPES [KAPPA_ATOMTYPES ...]] [-T TMIN] [-M TMAX]

Parse GROMACS processed topology file and perform desired scaling.

optional arguments:
  -h, --help            show this help message and exit
  -p TOPOL, --topol TOPOL
                        Processed Topology (default: processed.top)
  --show-molecules      Show molecule indices found in processed.top
  --show-molecule-atomtypes SHOW_MOLECULE_ATOMTYPES
                        Show molecule atomtypes found in processed.top
  --verbose
  -O OUTFILE, --outfile OUTFILE
                        Scaled topology base name (default: topol)
  -P OUTPUT, --output OUTPUT
                        Scaled topology output PATH (default: ./)
  -m METHOD, --method METHOD
                        Scaled topology method [ rest2 | ssrest3 ] (default: rest2)
  -H HOT_MOLECULES [HOT_MOLECULES ...], --hot-molecules HOT_MOLECULES [HOT_MOLECULES ...]
                        Space delimited entry of hot molecule selection (default: [0])
  -n NREPS, --nreps NREPS
                        Number of replicas (default: 10)
  -k KAPPA_MAX, --kappa-max KAPPA_MAX
                        Max Kappa Scaling (default: 1.06)
  --kappa-low-temp KAPPA_LOW_TEMP
                        Replicas at or above this temperature will have active solvent scaling (default:300)
  --kappa-atomtypes KAPPA_ATOMTYPES [KAPPA_ATOMTYPES ...]
                        Atomtypes to apply solvent scaling (default:['OW'])
  -T TMIN, --tmin TMIN  Base Temperature (default: 300)
  -M TMAX, --tmax TMAX  Max Temperature (default: 500)

Happy scaling!
```

### First query the molecules in the topology
```bash
username@computer] python3 src/repex_topology_parser.py -p tests/topology_files/test_topo/processed.top --show-molecules
0: Protein_chain_A
1: HxD
2: SOL
3: NA
4: CL
```

### Scaling with REST2 selecting 'Protein_chain_A' is simple:
```bash
username@computer] python3 src/repex_topology_parser.py -p tests/topology_files/test_topo/processed.top -p tests/topology_files/test_topo/processed.top -P tests/test_cli/ -m rest2 -H 0 -n 20 -O topol_rest2 --verbose
topology file read
Running rest2 scaling method
```

# For ssREST3, first to identify the atom(s) that will be scaled. For the a99disp-water model it is 
# important to note only the Oxygen heavy atom contains a non-zero LJ $\epsilon$ value and solvent scaling
# only applies to the LJ $\epsilon$ term. Therefore only the associated Oxygen atomtype is selected.
```bash
username@computer] python3 src/repex_topology_parser.py -p tests/topology_files/test_topo/processed.top --show-molecule-atomtypes 2
0: HW
1: MW
2: OW_tip4pd

# The command line ssREST3 scaling providing all options to cercumvent interactive mode:
```bash
username@computer] python3 src/repex_topology_parser.py -p tests/topology_files/test_topo/processed.top -p tests/topology_files/test_topo/processed.top -P tests/test_cli/ -m ssrest3 -H 0 -n 20 -k 1.1 --kappa-atomtypes OW_tip4pd -O topol_ssrest3 --verbose
topology file read
Running ssrest3 scaling method
```

### With interactive mode you can simply provide the topology, recognize the default values applied and 
### answer the questionaires to produce scaled topologies, however, this is discouraged for production runs. 
### Lastly, poor exchange ratios between specific replicas may occur and from within a jupyter notebook you 
### can supply the 'temps_opt':[ <temps(float)> ] to fine tune the temperature ladder. 
