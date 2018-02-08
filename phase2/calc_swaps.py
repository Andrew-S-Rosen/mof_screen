def update_calc(calc,calc_swaps):
#update calculator based on calc swaps
	for swap in calc_swaps:
		swap.replace(' ','')
		if swap == 'large_supercell':
			calc.special_params['lreal'] = 'Auto'
		elif swap == 'zbrent':
			calc.int_params['ibrion'] = 3
			calc.exp_params['ediff'] = 1e-6
			calc.int_params['nelmin'] = 8
			calc.int_params['iopt'] = 7
			calc.float_params['potim'] = 0
			calc.string_params['algo'] = 'Fast'
			calc_swaps.append('vtst')
		elif swap == 'dentet' or swap == 'grad_not_orth':
			calc.int_params['ismear'] = 0
			calc.string_params['algo'] = 'Fast'
		elif swap == 'edddav':
			calc.string_params['algo'] = 'All'
		elif swap == 'inv_rot_mat':
			calc.exp_params['symprec'] = 1e-8
		elif swap == 'subspacematrix' or swap == 'real_optlay' or swap == 'rspher' or swap == 'nicht_konv':
			calc.special_params['lreal'] = False
			calc.string_params['prec'] = 'Accurate'
		elif swap == 'tetirr' or swap == 'incorrect_shift':
			calc.input_params['gamma'] = True
		elif swap == 'rot_matrix':
			calc.input_params['gamma'] = True
			calc.int_params['isym'] = 0
		elif swap == 'pricel':
			calc.exp_params['symprec'] = 1e-8
			calc.int_params['isym'] = 0
		elif swap == 'amin':
			calc.float_params['amin'] = 0.01
		elif swap == 'pssyevx' or swap == 'eddrmm':
			calc.string_params['algo'] = 'Normal'
		elif swap == 'zheev':
			calc.string_params['algo'] = 'Exact'
		elif swap == 'elf_kpar':
			calc.int_params['kpar'] = 1
		elif swap == 'rhosyg':
			calc.exp_params['symprec'] = 1e-4
			calc.int_params['isym'] = 0
		elif swap == 'posmap':
			calc.exp_params['symprec'] = 1e-6
		elif 'sigma=' in swap:
			calc.float_params['sigma'] = float(swap.split('=')[-1])
		elif 'nbands=' in swap:
			calc.int_params['nbands'] = int(swap.split('=')[-1])
		elif 'potim=' in swap:
			calc.float_params['potim'] = float(swap.split('=')[-1])
		elif 'nsw=' in swap:
			calc.int_params['nsw'] = int(swap.split('=')[-1])
		elif 'nelm=' in swap:
			calc.int_params['nelm'] = int(swap.split('=')[-1])
		elif 'ibrion=' in swap:
			calc.int_params['ibrion'] = int(swap.split('=')[-1])
		elif 'istart=' in swap:
			calc.int_params['istart'] = int(swap.split('=')[-1])
		elif 'algo=' in swap:
			calc.string_params['algo'] = swap.split('=')[-1]
		elif 'iopt=' in swap:
			calc.int_params['iopt'] = int(swap.split('=')[-1])
		elif 'isif=' in swap:
			calc.int_params['isif'] = int(swap.split('=')[-1])
		elif 'lreal=' in swap:
			swap_val = swap.split('=')[1].lower()
			if swap_val == 'false':
				calc.special_params['lreal'] = False
			elif swap_val == 'auto':
				calc.special_params['lreal'] = 'Auto'
			elif swap_val == 'true':
				calc.special_params['lreal'] = True
	return calc