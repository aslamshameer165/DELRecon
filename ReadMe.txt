Main scripts/functions:
1. CT_MR_preprocess.m (Image preprocessing)
2. cluster_electrodes.m (Contact detection)
3. labeling.m (ITKSnap label file)

--- Run all the scripts in the above sequence to result in an electrode image with contact label file.
--- Download SPM12 toolbox (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and add it to the MATLAB path.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Input:
Path/<Subject_folder>/Imaging/

T1F.nii -- T1 MRI
CTpre.nii -- Pre-implantation CT
CTpost.nii -- Post-implantation CT
EleMap.txt -- Electrode and contact description (required for individual contact labeling)
%%%%%%%%%%%%%%%%%%%%%%%%
Output:

rmT1F.nii -- T1 MRI coregistered on CTpost.nii
rbrain_T1F_mask_ero.nii -- Smoothed brain mask
brain_masked_CTdiff.nii -- Brain masked CT difference image
initCL.mat -- Initial clustering results (Clustered electrode handle tip)
EleMap_NEW_updated.txt -- Updated EleMap.txt based on user defined values through GUI2
Contacts.mat -- Contact labeling information
all_electrodes.nii -- Electrode image (after contact detection)
Contact_labels.txt -- ITKSnap compatible contact label file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Dependencies for cluster_electrodes.m: 
- plotBrainandElec.m
- dataTipFormat.m
- clearCLUSTERs.m
- svdfit.m
- svdfitDist.m
- svd3D_1D.m
- svd1D_3d.m
- uptFigure.m
- uptFigure2.m
- uptFigure3.m
- uptFigure4.m
- clustering.m
- ndistance.m
- recluster.m
- inup.m
- indw.m
- yrange1.m
- yrange2.m
- adcontact.m
- rmcontact.m
- range1.m
- range2.m
- rst.m
