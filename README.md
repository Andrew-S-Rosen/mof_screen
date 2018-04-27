# PyMOFScreen
Python workflow for high-throughput DFT screening of MOFs

## Setup

1. PyMOFScreen requires Python 3.x (support for Python 2.x is not guaranteed). If you do not alerady have Python installed, the easiest option is to download the [Anaconda](https://www.anaconda.com/download/) distribution.

2. Install the most recent version of Pymatgen (e.g. via `pip install pymatgen`).

3. Install my custom branch of ASE (denoted rASE) at this link or patch `vasp.py` in `ase/ase/calculators/vasp/vasp.py` using the ``vasp.py` script located in the PyMOFScreen repository. The required modification adds support for checking if a VASP job has failed due to SCF convergence issues (via `atoms.calc.scf_converged`) or if it has reached the maximum number of ionic steps (via `atoms.calc.nsw_converged`).

4. VASP must be installed on your compute cluster and must be compiled with [VTSTools](http://theory.cm.utexas.edu/vtsttools/index.html).

5. Follow the `compute_environ` instructions below to set up PyMOFScreen to run on your compute environment.

## Ready-to-Run Examples

The main use of PyMOFScreen is to provide an easy-to-use DFT workflow built upon the Python packages ASE and Pymatgen for interfacing with VASP. Sample scripts are provided in `/examples` that can be used to: 1) perform a full volume relaxation on a series of CIFs (`volume_relaxation.py`); 2) perform an ionic relaxation on a series of CIFs (`ionic_relaxation.py`). 

## screener

The main tool to initialize a screening workflow is the `pymofscreen.screener` clas, which is described below. At the bare minimum, must provide it the path to where the CIF files are located (`mofpath`), the base directory where the DFT screening results should be stored (`basepath`), and the name of the job submission script (`submit_script`). All the re and all errors to `basepath/errors`. sults will be written to the folder `basepath/results` and all errors to `basepath/errors`. By default, `screener` assumes that you want the CIF file to be converted to the Niggli-reduced unit cell, which can be changed via the `niggli` keyword. It also assumes that the stdout file is simply the basename of the Python submission script with a `.out` extension, which can be changed via `stdout_file`. Finally, there are the k-point details. By default, the workflow assumes you want a low-accuracy k-point per atom (KPPA) of 100 and a high-accuracy KPPA of 1000, which can be changed by changing `kppas`. Generally, `kpts_path` does not need to be altered unless you wish to manually specify the k-points for each CIF.

```python
class screener():
	"""
	This class constructs a high-throughput screening workflow
	"""
	def __init__(self,mofpath,basepath,kpts_path='Auto',kppas=None,niggli=True,
		submit_script='sub_screen.job',stdout_file=None):
		"""
		Initialize variables that should be used on all MOFs in a database
		Essential args:
			mofpath (string): path to the directory containing the CIF files
      
			basepath (string): path to the base directory for the DFT screening
			
      kppas (list of ints): the low-accuracy and high-accuracy KPPAs to use if
      kpts_path == 'Auto' (defaults to [100, 1000])
			
      niggli (bool): True if Niggli-reduction should be performed
			
      submit_script (string): name of job submission script
      (defaults to sub_screen.job)
			
      stdout_file (string): name of the stdout file (defaults to the
      basename of the .py script with a .out extension)
		"""
```

Within the `screener` class is a function named `run_screen`, which is described below. It informs the `screener` what type of job should be run and on what CIF file. Generally, only two parameters ever need to be changed: the name of the CIF file (`cif_file`) and the type of job to be run (`mode`), which can be either `volume` or `ionic`.`

```python
	def run_screen(self,cif_file,mode,spin_levels=None,acc_levels=None,calcs=calcs):
		"""
		Run high-throughput ionic relaxations
		Args:
			cif_file (string): name of CIF file
		Returns:
			mofs (list of ASE Atoms objects): ASE Atoms objects for optimized
			MOF given by cif_file for each spin_level
		"""
```
## Example

Now, that was a bit abstract. It's pretty easy in practice though! The `volume_relaxation.py` example is reproduced below. Note that the `screener` is first initialized as previously described, and then every CIF file is looped over and `run_screen` is performed. Note that there is a function `pymofscreen.cif_handler.get_cif_files`, which will automatically make a list of the names of all CIF files in `mofpath`.

```python
from pymofscreen.cif_handler import get_cif_files
from pymofscreen.screen import screener

#Run screening analysis
mofpath = '/projects/p30148/vasp_jobs/structures/CoRE1-DFT-OMS/'
basepath = '/projects/p30148/vasp_jobs/MOFs/testing/'
submit_script = 'sub_screen.job'
cif_files = get_cif_files(mofpath)
s = screener(mofpath,basepath,submit_script=submit_script)
for cif_file in cif_files:
	mof = s.run_screen(cif_file,'volume')
```
## Compute Environments

Every compute environment is different, with different ways to run VASP and different job submission systems. To address this, the `pymofscreen.compute_environ` file must be altered before the first use. In `compute_environ.get_nprocs`, the variabels `nprocs` and `ppn` must be integers representing total number of processors and the processors per node, respectively. This is set up to correctly parse the information from a MOAB-style submit script by default (this is where the variable `submit_script` comes into play) but must be changed for your system. In `compute_environ.choose_vasp_version`, you must inform PyMOFScreen how to properly run VASP on the given machine (i.e. how to correctly set up ASE's `run_vasp.py` file). See [here](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html) for details of the `run_vasp.py` file that ASE required. Generally, you will just need to tell `choose_vasp_version` the module name of VASP, the names of the executables for the gamma-point and standard versions of the VTST-enabled VASP builds, and how to run VASP on the given machine.

To help you get started, for both functions, examples are provided for Quest at Northwestern, Cori at NERSC, Stampede2 at TACC, and Thunder at AFRL. 
