import os
from pymofscreen.writers import pprint
from pymofscreen.kpts_handler import get_kpts
from pymofscreen.calculators import defaults
from pymofscreen.screen_phases import workflows
from pymofscreen.janitor import prep_paths

class screener():
	"""
	This class constructs a high-throughput screening workflow
	"""
	def __init__(self,mofpath,basepath,ads_species=None,kpts_path=None,kpts_method=None,niggli=False,
		spin_levels=['spin1','spin2'],acc_levels=['scf_test','isif2_lowacc',
		'isif2_medacc','final','final_spe'],submit_script='sub_screen.job',stdout_file='driver.out'):
		self.mofpath = mofpath
		self.basepath = basepath
		self.ads_species = ads_species
		self.submit_script = submit_script
		self.stdout_file = stdout_file
		self.kpts_path = kpts_path
		self.kpts_method = kpts_method
		self.niggli = niggli
		self.spin_levels = spin_levels
		self.acc_levels = acc_levels
		if self.kpts_path == self.kpts_method:
			raise ValueError('Set one of kpts or kpts_method')
		prep_paths(basepath)

	def run_ads_screen(self,cif_file):

		#Run high-throughput screening
		basepath = self.basepath
		spin_levels = self.spin_levels

		#Make sure MOF isn't running on other process
		working_cif_path = os.path.join(basepath,'working',cif_file)
		if os.path.isfile(working_cif_path) == True:
			pprint('SKIPPED: Running on another process')
			return None

		#Get the kpoints
		kpts_lo, gamma = get_kpts(self,cif_file,defaults['kppa_lo'])
		kpts_hi, gamma = get_kpts(self,cif_file,defaults['kppa_hi'])
		defaults['gamma'] = gamma
		defaults['kpts_lo'] = kpts_lo
		defaults['kpts_hi'] = kpts_hi

		mofs = []
		#for each spin level, optimize the structure
		for i, spin_level in enumerate(spin_levels):

			#***********PREP FOR RUN***********
			if spin_level != spin_levels[0]:
				prior_spin = spin_levels[i-1]
			else:
				prior_spin = None
			wf = workflows(self,defaults,cif_file,spin_level,prior_spin)

			#***********SCF TEST************
			mof = wf.scf_test()
			if mof is None:
				return None

			#***********ISIF 2 (lowacc)************
			mof = wf.isif2_lowacc()
			if mof is None:
				return None

			#***********ISIF 2 (medacc)************
			mof = wf.isif2_medacc(mof)
			if mof is None:
				return None

			#***********ISIF 2 (final)************
			mof = wf.isif2_final(mof)
			if mof is None:
				return None
				
			#***********Final SPE************
			mof, skip_spin2 = wf.final_spe(mof)
			if mof is None:
				return None

			#***********SAVE and CONTINUE***********
			mofs.append(mof)
			if skip_spin2 == True:
				pprint('Skipping '+spin_levels[i+1]+' run')
				return mofs

		return mofs