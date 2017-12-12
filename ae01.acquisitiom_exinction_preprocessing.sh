#!/bin/tcsh
set pdir=/media/Lapis/DougDiss_MVCS/ACQ
set qdir=/media/Lapis/DougDiss_MVCS/EXT
source /usr/share/Modules/init/csh
module load afni_6.9.15

#Run ACQ & EXT preprocessing for all subjects
cd ${pdir}
#	foreach direct ( "ACQ" "EXT" )
	foreach subject (*)
		set subdir=${pdir}/${subject} ;
		set subdirext=${qdir}/${subject} ;
		cd ${subdir} ;
		echo "Preprocessing RS for subject ${subject}" ;
	

##########################################################################################
##############################	LET'S START	##########################################
##########################################################################################

	#!/bin/tcsh -xef
	#make sure that you're in tcsh
	#make sure you're in the data parent directory

# script setup

#  take note of the AFNI version
afni -ver

#  the user may specify a single subject to run with
if ( $#argv > 0 ) then
    set subj = $argv[1]
else
    set subj = ${subject}
endif

# assign output directory name
set output_dir = $subdir/$subj.ACQ.results
set output_dir_ext = $subdirext/$subj.EXT.results

# verify that the results directory does not yet exist
if ( -d $output_dir ) then
    echo output dir "$subj.ACQ.results" already exists
    rm -r $output_dir
endif

if ( -d $output_dir_ext ) then
    echo output dir "$subj.EXT.results" already exists
    rm -r $output_dir_ext
endif

#  set list of runs - onl 1 for both ACQ & EXT
set runs = (`count -digits 2 1 1`)
# create results and stimuli directories
mkdir $output_dir
mkdir $output_dir/stimuli

mkdir $output_dir_ext
mkdir $output_dir_ext/stimuli

# copy anatomy to results dir
3dcopy ${subdir}/1${subj}.anat.al+orig \
    $output_dir/1${subj}.anat.al

3dcopy ${subdirext}/1${subj}.anat.al+orig \
    $output_dir_ext/1${subj}.anat.al

# ============================ tcat ============================
# apply 3dTcat to copy input dsets to results dir, while
# removing the first 4 TRs

3dTcat -prefix $output_dir/pb00.$subj.r01.tcat                              \
    ${subdir}/1F1${subj}3d+orig'[4..$]'
3dTcat -prefix $output_dir_ext/pb00.$subj.r01.tcat                              \
    ${subdirext}/2F3${subj}3d+orig'[4..$]'

# and make note of repetitions (TRs) per run
set tr_counts = ( 286 )

# -------------------------------------------------------
# enter the results directory (can begin processing data)
cd $output_dir

# ========================== auto block: outcount ==========================
# data check: compute outlier fraction for each volume
touch out.pre_ss_warn.txt
foreach run ( $runs )
    3dToutcount -automask -fraction -polort 4 -legendre                     \
                pb00.$subj.r$run.tcat+orig > outcount.r$run.1D
	3dToutcount -automask -fraction -polort 4 -legendre                     \
                $output_dir_ext/pb00.$subj.r$run.tcat+orig > $output_dir_ext/outcount.r$run.1D

    # outliers at TR 0 might suggest pre-steady state TRs
    if ( `1deval -a outcount.r$run.1D"{0}" -expr "step(a-0.4)"` ) then
        echo "** TR #0 outliers: possible pre-steady state TRs in run $run" \
            >> out.pre_ss_warn.txt
    endif
end

# catenate outlier counts into a single time series
cat outcount.r*.1D > outcount_rall.1D
cat $output_dir_ext/outcount.r*.1D > $output_dir_ext/outcount_rall.1D

# ================================= tshift =================================
# time shift data so all slice timing is the same 
foreach run ( $runs )
	echo "** Time shifting subject ${subj} run ${run}"
    3dTshift -tzero 0 -quintic -prefix pb03.$subj.r$run.tshift \
             pb00.$subj.r01.tcat+orig
    3dTshift -tzero 0 -quintic -prefix $output_dir_ext/pb03.$subj.r$run.tshift \
             $output_dir_ext/pb00.$subj.r01.tcat+orig
end
	
cd - 
cd ${output_dir}
#This is the actual command which is much more compact than individual processing blocks. Make sure to turn on volreg option since you have not motion corrected yet.

#I've decided to run alignment sepaately for AXQ and EXT because EXT needs its own motion correction which it does not get as a child EPI.
	echo "** Motion Correcting subject ${subj} ACQ and aligning all runs to run 1 which is aligned to anatomical 1"
	align_epi_anat.py \
		-overwrite \
		-anat 1${subj}.anat.al+orig \
		-epi pb03.${subj}.r01.tshift+orig  \
		-epi_base 0 \
		-volreg_method 3dAllineate \
		-master_epi pb03.${subj}.r01.tshift+orig \
		-tshift off \
		-deoblique on \
		-epi2anat \
		-giant_move \
		-save_vr \
		-save_skullstrip \
		-suffix _ALIGN

#EXT alignment
cd ${output_dir_ext}
	echo "** Motion Correcting subject ${subj} EXT and aligning all runs to run 1 which is aligned to anatomical 1"
	align_epi_anat.py \
		-overwrite \
		-anat 1${subj}.anat.al+orig \
		-epi pb03.${subj}.r01.tshift+orig  \
		-epi_base 0 \
		-volreg_method 3dAllineate \
		-tshift off \
		-deoblique on \
		-epi2anat \
		-giant_move \
		-resample off \
		-save_vr \
		-save_skullstrip \
		-suffix _ALIGN

cd ${output_dir}
cat pb03.$subj.r*.tshift_vr_motion.1D > dfile_rall.1D
cat $output_dir_ext/pb03.$subj.r*.tshift_vr_motion.1D > $output_dir_ext/dfile_rall.1D


1d_tool.py -infile dfile_rall.1D  -set_run_lengths 286 \
           -derivative  -collapse_cols euclidean_norm               \
           -write motion_${subj}_enorm.1D
1d_tool.py -infile $output_dir_ext/dfile_rall.1D  -set_run_lengths 286 \
           -derivative  -collapse_cols euclidean_norm               \
           -write $output_dir_ext/motion_${subj}_enorm.1D

########REST OF THIS IS OPTIONAL######## 
##Can use pb03.$subj.r05.tshift_ALIGN files if rest of processing looks weird##

# ================================== mask ==================================
# Output from previous step is pb03.$subj.r$run.tshift_ALIGN+orig\
# create 'full_mask' dataset (union mask)
foreach run ( $runs )
echo "** Automasking subject ${subj} run ${run} "
    3dAutomask -dilate 1 -prefix rm.mask_r$run pb03.$subj.r$run.tshift_ALIGN+orig
    3dAutomask -dilate 1 -prefix $output_dir_ext/rm.mask_r$run $output_dir_ext/pb03.$subj.r$run.tshift_ALIGN+orig
end

# get mean and compare it to 0 for taking 'union'
3dMean -datum short -prefix rm.mean rm.mask*.HEAD
3dcalc -a rm.mean+orig -expr 'ispositive(a-0)' -prefix full_mask.$subj
3dMean -datum short -prefix $output_dir_ext/rm.mean rm.mask*.HEAD
3dcalc -a $output_dir_ext/rm.mean+orig -expr 'ispositive(a-0)' -prefix $output_dir_ext/full_mask.$subj

# ---- create subject anatomy mask, mask_anat.$subj+orig ----
#      (resampled from aligned anat which is just original anat because we did an EPI -> ANAT alignment)
3dresample -master full_mask.$subj+orig -input 1${subj}.anat.al_ns+orig \
           -prefix rm.resam.anat
3dresample -master $output_dir_ext/full_mask.$subj+orig -input $output_dir_ext/1${subj}.anat.al_ns+orig \
           -prefix $output_dir_ext/rm.resam.anat

# convert to binary anat mask; fill gaps and holes
3dmask_tool -dilate_input 5 -5 -fill_holes -input rm.resam.anat+orig  \
            -prefix mask_anat.$subj
3dmask_tool -dilate_input 5 -5 -fill_holes -input $output_dir_ext/rm.resam.anat+orig  \
            -prefix $output_dir_ext/mask_anat.$subj

# compute overlaps between anat and EPI masks
3dABoverlap -no_automask full_mask.$subj+orig mask_anat.$subj+orig    \
            |& tee out.mask_overlap.txt
3dABoverlap -no_automask $output_dir_ext/full_mask.$subj+orig $output_dir_ext/mask_anat.$subj+orig    \
            |& tee $output_dir_ext/out.mask_overlap.txt

# ---- segment anatomy into classes CSF/GM/WM ----
# Create an 'anat_final' dataset which is he anat everything is aligned to
3dcopy 1$subj.anat.al_ns+orig anat_final.$subj
3dcopy $output_dir_ext/1$subj.anat.al_ns+orig $output_dir_ext/anat_final.$subj

end
