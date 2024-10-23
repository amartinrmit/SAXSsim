import subprocess

#
# A script to choose parameters and run saxs.py
#

# the .pdb file with the atomic structure
f = 'TBAB500EtOH1000_1g_100ns_ff_noH.pdb'

# text to append to every output file (change this to avoid overwriting previous calculations)
suffix = "_norf_TEST"

# the location where the output will be be saved
outpath = path = f"/Users/andrewmartin/cloudstor/Work/Research/SF/Projects/2024/dynasor/output2/"

# number of samples used to estimate the background from the mean density in the volume
saxs_box_nx = 500

# select a ragne for r-values to normalise the background before subtraction
saxs_norm_rmin, saxs_norm_rmax = 5, 10

# run saxs.py
tag = f[:-4]+suffix
flags = f'--pdbname {f} --tag {tag} --outpath {outpath} --saxs_box_nx {saxs_box_nx} --saxs_norm_rmin {saxs_norm_rmin} --saxs_norm_rmax {saxs_norm_rmax}'
result = subprocess.run('python /Users/andrewmartin/cloudstor/Work/Research/SF/codes/code_projects/saxs_tmp/saxs.py -c config_saxs.txt '+flags,shell=True)
