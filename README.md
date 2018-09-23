# PyMOFScreen
Python workflow for high-throughput DFT screening of MOFs using VASP. Relevant details for the code can be found in the following manuscript (once published):

A.S. Rosen, J.M. Notestein, R.Q. Snurr. "Identifying Promising Metal-Organic Frameworks for Heterogeneous Catalysis via High-Throughput Periodic Density Functional Theory." In preparation. 

## What is PyMOFScreen?

High-throughput computational catalysis involving MOFs is a tricky business. Their large unit cells, diverse structures, and widely varying compositions make it challenging to achieve both a robust and high-performing workflow with little human interactions. PyMOFScreen solves this problem through multi-stage structural optimizations, a robust selection of optimization algorithms that are chosen on-the-fly, automatic error-handling, and much more. In the Snurr group, we have used PyMOFScreen to screen hundreds of MOFs using periodic DFT in a fully automated fashion.

## Ready-to-Run Examples

To get started, sample scripts are provided in `/examples` and include:
1. `volume_relaxation.py`. Perform a full volume relaxation on a database of MOF CIFs.
2. `ionic_relaxation.py`. Perform an ionic relaxation on a database of MOF CIFs.

## The screener

The main tool to initialize a screening workflow is the `pymofscreen.screener` class, which is described below. At the bare minimum, you must provide it the base directory where the DFT screening results should be stored (`basepath`), the path to where the MOF CIF files are located (`mofpath`), the name of the job submission script (`submit_script`), and the name of the standard output (`stdout_file`). All the results are stored in `basepath/results`, and any errors are stored in `basepath/errors`. The `screener` also requires that the user specify the k-points per atom (KPPA) that should be used. By default, it will use 100 KPPA and 1000 KPPA for the low- and high-accuracy phases of the workflow, respectively. Generally, `kpts_path` does not need to be altered unless you wish to manually specify the k-points for each CIF.

```python
class screener():
	"""
	This class constructs a high-throughput screening workflow
	"""
	def __init__(self,basepath,mofpath=None,kpts_path='Auto',kppas=None,
		submit_script=None,stdout_file=None):
		"""
		Initialize variables that should be used on all MOFs in a database
		Args:
			basepath (string): path to the base directory for the DFT screening

			mofpath (string): path to the directory containing the CIF files

			kpts_path (string): can be either 'Auto' for an automatic generation
			of the kpoints based on KPPAs or a string representing the path to a
			text file with all the kpoint information (refer to examples/kpts.txt)

			kppas (list of ints): KPPAs to use if kpts_path == 'Auto' (defaults
			to kppas = [100, 1000] for 100 and 1000 KPPA for the low and high
			accuracy runs)

			submit_script (string): path to job submission script

			stdout_file (string): path to the stdout file (defualts to the 
			name of the Python job with a .out extension instead of .py)
		"""
```

Within the `screener` class is a function named `run_screen`, which is described below. It informs the `screener` what type of job should be run and on what CIF file. Generally, two parameters need to be changed: the name of the CIF file (`cif_file`) and the type of job to be run (`mode`), which can be either `volume` or `ionic`. By default `spin_levels` parameter is set to `[spin1,spin2]` such that a high-spin and then low-spin job is performed. By default, `acc_levels` is set to `['scf_test','isif2_lowacc','isif2_medacc','isif2_highacc','final_spe']` such that this sequence of jobs is performed (as discussed below). Also, `niggli` specifies whether the unit cell should be Niggli-reduced and defaults to `True`.

```python
def run_screen(self,cif_file,mode,spin_levels=None,acc_levels=None,niggli=True,calcs=calcs):
	"""
	Run high-throughput ionic or volume relaxations
	Args:
		cif_file (string): name of CIF file

		mode (string): 'ionic' or 'volume'

		spin_levels (list of strings): spin states to consider (defaults
		to ['spin1','spin2'])

		acc_levels (list of strings): accuracy levels to consider (defaults
		to ['scf_test','isif2_lowacc','isif2_medacc','isif2_highacc','final_spe'])

		niggli (bool): True/False if Niggli-reduction should be done (defaults
		to niggli=True)

		calcs (function): function to call respective calculator (defaults to
		automatically importing from pymofscreen.default_calculators.calcs)

	Returns:
		best_mof (ASE Atoms objects): ASE Atoms object for optimized MOF
	"""
```
## Example

Now, that was a bit abstract. It's pretty easy in practice though! A minimal example for performing a volume relaxation is shown below. There is a function `pymofscreen.cif_handler.get_cif_files`, which will automatically make a list of the names of all CIF files in `mofpath`. Then, the `screener` is first initialized, and `run_screen` is performed for every CIF file. 

```python
from pymofscreen.cif_handler import get_cif_files
from pymofscreen.screen import screener
#Set up paths
mofpath = 'PathToCIFs'
basepath = 'PathToStoreResults'
submit_script = 'PathToSubmitScript'
#Get CIF files
cif_files = get_cif_files(mofpath)
#Construct screener object
s = screener(basepath,mofpath,submit_script=submit_script)
#Run screening
for cif_file in cif_files:
	mof = s.run_screen(cif_file,'volume')
```

## Defaults

Of course, it is essential to specify default parameters that should be used in VASP, such as the exchange-correlation functional, convergence criteria, and so on. This is done by importing `pymofscreen.default_calculators.defaults` and making modifications to the default parameters in the `defaults` dictionary. An example is shown below.

```python
from pymofscreen.cif_handler import get_cif_files
from pymofscreen.screen import screener
from pymofscreen.default_calculators import defaults
#Set up paths
mofpath = 'PathToCIFs'
basepath = 'PathToStoreResults'
submit_script = 'PathToSubmitScript'
#Define defaults
defaults['xc'] = 'M06L'
defaults['ivdw'] = 11
defaults['ediffg'] = -0.02 #and so on...
#Get CIF files
cif_files = get_cif_files(mofpath)
#Construct screener object
s = screener(basepath,mofpath,submit_script=submit_script)
#Run screening
for cif_file in cif_files:
	mof = s.run_screen(cif_file,'volume')
```

The parameters in the `defaults` dictionary are used in the `pymofscreen.default_calculators.calcs` function, which we suggest looking at before running PyMOFScreen for the first time. The `pymofscreen.default_calculators.calcs` function defines each job type previously specified in `acc_levels` within `run_screen`. For instance, it defines `isif2_lowacc` as an low accuracy ionic relaxation and `final_spe` as a high accuracy, single point energy calculation using the parameters stored in `defaults`. The job specifications and parameters can be freely changed using any of [ASE's parameters for VASP](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html).

## Setup

### Installing PyMOFScreen

1. PyMOFScreen requires [Python](https://www.python.org/) 3.6 or newer. If you do not already have Python installed, the easiest option is to download the [Anaconda](https://www.anaconda.com/download/) distribution.
2. Download or clone the PyMOFScreen repository and run `pip install -r requirements.txt` followed by `pip install .` from the PyMOFScreen base directory. This will install PyMOFScreen and the required dependencies (rASE, Pymatgen).

### Required Dependencies

PyMOFScreen requires the following Python packages. Both are installed by using `pip install -r requirements.txt`.
1. [Pymatgen](http://pymatgen.org/) 2018.5.22 or newer.
2. A slightly modified build of [ASE](https://wiki.fysik.dtu.dk/ase/) 3.16.2 or newer. The required modification adds support for checking if a VASP job has failed due to SCF convergence issues (via `atoms.calc.scf_converged`) and if it has reached the maximum number of ionic steps (via `atoms.calc.nsw_converged`). The custom build, denoted rASE, can be found [at this link](https://github.com/arosen93/rASE). Alternatively, you can directly patch `vasp.py` in `ase/ase/calculators/vasp/vasp.py` using the `vasp.py` script located in the PyMOFScreen `patches` directory. Regardless, ensure that the `VASP_PP_PATH` environment variable is set according to the details [here](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html).

PyMOFScreen also requires that VASP is installed on your compute cluster. The VASP build must be compiled with [VTSTools](http://theory.cm.utexas.edu/vtsttools/index.html) and must include both gamma-point only and standard builds. Follow the `compute_environ` instructions below to set up PyMOFScreen to run on your compute environment.

### Compute Environments

Every compute environment is unique, with different ways to run VASP and different job submission systems. To address this, the `pymofscreen/compute_environ` file must be altered. Templates have been provided for Quest at Northwestern, Cori at NERSC, Stampede2 at TACC, and Thunder at AFRL that you can uncomment as needed. If you are using another machine, follow the instructions below.

In `compute_environ.get_nprocs`, the variable `nprocs` must be able to determine the number of processors for a given job from the submission script. In `compute_environ.choose_vasp_version`, you must inform PyMOFScreen how to properly run VASP on the given machine (i.e. how to correctly set up ASE's `run_vasp.py` file). See [here](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html) for details of the `run_vasp.py` file that ASE requires. Generally, you will just need to tell `choose_vasp_version` the module name of VASP, the names of the executables for the gamma-point and standard versions of the VTST-enabled VASP builds, and how to run VASP on the given machine.