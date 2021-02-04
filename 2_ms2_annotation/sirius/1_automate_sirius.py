#!/usr/bin/python

"""
Usage: $ ./automate_sirius.py

Goal: Automate Sirius annotation of ms2 spectra

---
NGG.
"""

import glob,os,re,sys,shutil,subprocess

outdir = 'SIRIUS_output'
indir = 'groups'

files = glob.glob(indir+'/*.mgf') 

for file in files:
    group = file.split('/')[-1].split('.')[0]
    print('###### Working with {}'.format(group))
    print()
    subprocess.call(['./5_sirius_cli_run.sh', os.path.join(outdir, group), file])
