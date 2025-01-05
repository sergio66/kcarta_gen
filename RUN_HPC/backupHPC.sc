## /bin/cp -a *.m *.nml *.sbatch *.sc *eadme* .

/bin/cp -a ../WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/*.m      .
/bin/cp -a ../WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/*.nml    .
/bin/cp -a ../WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/*.sbatch .
/bin/cp -a ../WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/*.sc     .
/bin/cp -a ../WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/*eadme*  .

## /bin/cp -aR ../WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/*/*.m  OUTPUT_mfiles/.

#rsync -am --include='*.m' --include='*/' --exclude='*' testdir/ testdir2/
#where testdir is the top level of the tree you want to copy *.m files
#from and testdir2 is a directory you mkdir to hold the results 10:38

## first do dry run
rsync -avm --dry-run --include='*.m' --include='*/' --exclude='*' ../WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ OUTPUT_mfiles/
## then actual run
rsync -avm           --include='*.m' --include='*/' --exclude='*' ../WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ OUTPUT_mfiles/
