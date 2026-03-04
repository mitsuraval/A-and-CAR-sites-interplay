# This script writes the cpptraj input scripts for out hydrogen bond analysis

#################
# JOB VARIABLES #
#################

start_expt = 1
end_expt = 30

data_dir = '/home66/mraval/tRNAmod/correlations/AA_G34/AACGCU/NEUTRAL/NEUTRAL_'
# use solvated (wat) prmtop
prmtop_path = '/home66/mraval/tRNAmod/34GUU36/AA_G34/AACGCU/TLEAP/5JUP_N2_wat.prmtop'
trajin_path = '/home66/mraval/tRNAmod/34GUU36/AA_G34/AACGCU/NEUTRAL/NEUTRAL_'

# ns to remove from analysis based on RMSD
# usually 20ns
cutoff = 20

frames_per_ns = 100

####################
# CALCULATE CUTOFF #
####################

# Frame where the simulation stabilizes
stable_cutoff = str(frames_per_ns*cutoff)

###############
# WRITE FILES #
###############

# files to submit to scheduler
outfile = open('run_cpptraj_3a_avgHbond.sh','w')
outfile.write('#!/bin/bash\n' +\
    '#SBATCH --job-name="3a_avgHbond_CPU"\n' +\
    '#SBATCH --output=out_3a_avgHbond\n' +\
    '#SBATCH --error=err_3a_avgHbond\n' +\
    '#SBATCH --ntasks=1\n' +\
	'#SBATCH --cpus-per-task=4\n' +\
    '#SBATCH -N 1\n' +\
    '#SBATCH --partition=exx96\n' +\
    '#SBATCH -B 1:1:1\n' +\
    
    '# env\n' +\
    'export PATH=/share/apps/CENTOS7/gcc/6.5.0/bin:$PATH\n' +\
    'export LD_LIBRARY_PATH=/share/apps/CENTOS7/gcc/6.5.0/lib64:$LD_LIBRARY_PATH\n\n' +\

    '### AMBER22 ###\n' +\
    'source /share/apps/CENTOS7/amber/amber22/amber.sh\n\n' +\
    
    '### JOB-SPECIFIC COMMANDS ###\n\n' +\
    
	'# call cpptraj in Amber22\n' +\
	'for i in {' + str(start_expt) + '..' + str(end_expt) + '}; do\n' +\
	'	cd ' + data_dir + '$i/\n' +\
	'	cpptraj -i cpptraj_3a_avgHbond.in\n' +\
	'done\n\n')
outfile.close()

for i in range(start_expt,end_expt+1):
	outfile = open(data_dir + str(i) + '/cpptraj_3a_avgHbond.in', 'w')
	outfile.write('\n' +\
		'parm ' + prmtop_path + ' [modi3]\n' +\
		'trajin ' + trajin_path + str(i) + '/mdcrd_nd_' + str(i) + ' ' + ' parm [modi3]\n' +\
		'autoimage\n\n' +\
		
		'### In-Registration ###\n' +\
		'# Anticodon position 1/codon position 1\n' +\
		'hbond nhb_AVE_410_all_423_all :410|:423 out nhb_AVE_410_all_423_all.dat\n' +\
		'# Anticodon position 2/codon position 2\n' +\
		'hbond nhb_AVE_409_all_424_all :409|:424 out nhb_AVE_409_all_424_all.dat\n' +\
		'# Anticodon position 3/codon position 3\n' +\
		'hbond nhb_AVE_408_all_425_all :408|:425 out nhb_AVE_408_all_425_all.dat\n' +\
		'# C1054 to +1 N1\n' +\
		'hbond nhb_AVE_94_all_426_all :94|:426 out nhb_AVE_94_all_426_all.dat\n' +\
		'# A1196 to +1 N2\n' +\
		'hbond nhb_AVE_127_all_427_all :127|:427 out nhb_AVE_127_all_427_all.dat\n' +\
		'# R146 to +1 N2\n' +\
		'hbond nhb_AVE_240_all_427_all :240|:427 out nhb_AVE_240_all_427_all.dat\n\n' +\
		
		'### Cross-Registration ###\n' +\
		'# Anticodon position 1/codon position 2\n' +\
		'hbond nhb_AVE_410_all_424_all :410|:424 out nhb_AVE_410_all_424_all.dat\n' +\
		'# Anticodon position 2/codon position 1\n' +\
		'hbond nhb_AVE_409_all_423_all :409|:423 out nhb_AVE_409_all_423_all.dat\n' +\
		'# Anticodon position 2/codon position 3\n' +\
		'hbond nhb_AVE_409_all_425_all :409|:425 out nhb_AVE_409_all_425_all.dat\n' +\
		'# Anticodon position 3/codon position 2\n' +\
		'hbond nhb_AVE_408_all_424_all :408|:424 out nhb_AVE_408_all_424_all.dat\n' +\
		'# Anticodon position 3 to N1\n' +\
		'hbond nhb_AVE_408_all_426_all :408|:426 out nhb_AVE_408_all_426_all.dat\n' +\
		'# C1054 to codon position 3\n' +\
		'hbond nhb_AVE_94_all_425_all :94|:425 out nhb_AVE_94_all_425_all.dat\n' +\
		'# C1054 to +1 N2\n' +\
		'hbond nhb_AVE_94_all_427_all :94|:427 out nhb_AVE_94_all_427_all.dat\n' +\
		'# C1054 to +1 N3\n' +\
		'hbond nhb_AVE_94_all_428_all :94|:428 out nhb_AVE_94_all_428_all.dat\n' +\
		'# A1196 to +1 N1\n' +\
		'hbond nhb_AVE_127_all_426_all :127|:426 out nhb_AVE_127_all_426_all.dat\n' +\
		'# A1196 to +1 N3\n' +\
		'hbond nhb_AVE_127_all_428_all :127|:428 out nhb_AVE_127_all_428_all.dat\n' +\
		'# R146 to +1 N1\n' +\
		'hbond nhb_AVE_240_all_426_all :240|:426 out nhb_AVE_240_all_426_all.dat\n' +\
		'# R146 to +1 N3\n' +\
		'hbond nhb_AVE_240_all_428_all :240|:428 out nhb_AVE_240_all_428_all.dat\n\n')
	
	outfile.close()
