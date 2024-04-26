%% Preprocessing of CT and MRI images in preparation for DELRecon.
% This script implements the Module 1 of DELRecon
% Script assumes that scans are named T1F (T1 MRI), CTpre (pre-implantation CT) 
% and CTpost (post-implantation CT) are in NIfTI format in
% subdirectory /path/<subject_folder>/Imaging/
% Preprocessing steps include: Co-registration of pre-surgical CT and T1 MRI to post-surgical CT, tissue segmentation of MRI
% and creation of brain-masked CT difference image (CTdiff).
% Download and add path to the spm12


clear all
close all
%%  Define path to the data
basedir = 'E:\PROJECT\SEEG segmentation\Edited\DELRecon\'; % Path to the subject folder
sub = 'subject'; % subject folder to be processed
cd(basedir);
subjdir = strcat(basedir,sub,'\Imaging\'); %path to the data to be processed
cd(subjdir)
 

%%      Coregister CTpre to CTpost (using SPM12)  

% Define the coregistration parameters
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[subjdir 'CTpost.nii,1']}; %post-implantation CT
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[subjdir 'CTpre.nii,1']}; %pre-implantation CT
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r'; % save the coregistered 'CTpre' as 'rCTpre' 

% Run the coregistration
spm_jobman('run', matlabbatch);
clear matlabbatch 

%% Segment the T1 MRI into different tissue types

% Define the segmentation parameters
matlabbatch{1}.spm.tools.preproc8.channel.vols = {[subjdir 'T1F.nii,1']}; % pre-implantation T1 MRI
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 1];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spm('Dir') '\tpm\TPM.nii,1']}; % define tissue class for each tissue type through TPM.nii in SPM folder
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 1;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0]; % save only native tissue type
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spm('Dir') '\tpm\TPM.nii,2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 1;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spm('Dir') '\tpm\TPM.nii,3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spm('Dir') '\tpm\TPM.nii,4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[spm('Dir') '\tpm\TPM.nii,5']};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[spm('Dir') '\tpm\TPM.nii,6']};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 1;
matlabbatch{1}.spm.tools.preproc8.warp.cleanup = 1;
matlabbatch{1}.spm.tools.preproc8.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.fwhm = 0;
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1];

% Run the segmentation
spm_jobman('run', matlabbatch);
clear matlabbatch

%% Create the brain mask with grey matter, white matter and CSF

clear Vi Vo c V
Vi = strvcat([subjdir 'c1T1F.nii'], [subjdir 'c2T1F.nii'], [subjdir 'c3T1F.nii']); % grey matter (c1T1F.nii), white matter (c2T1F.nii) and CSF (c3T1F.nii) tissue maps for brain mask
Vo = [subjdir 'brain_T1F_mask.nii']; % The brain mask 
spm_imcalc(Vi, Vo, 'i1+i2+i3'); % Creating the brain mask by adding grey matter, white matter and CSF

% Dilate and erode the brain mask using spherical neighborhood to fill in
% all the voids in the brain mask which may be due to tissue
% misclassification from segmentation
V=spm_vol(Vo);
c=spm_read_vols(V); % Load the brain mask
[xx,yy,zz] = ndgrid(-5:5);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 4.0; % define a three-dimensional spherical neighborhood of 2 unit radius
c=imdilate(c,nhood);
c=imdilate(c,nhood);
c = imerode(c,nhood);
c = imerode(c,nhood)>0.9;
% c = imerode(c,nhood)>0.9;
V.fname='brain_T1F_mask_ero.nii';
spm_write_vol(V,c); % save the smoothed brain mask with CSF
clear br c1 c2 c3 c V Vi Vo matlabbatch

%% Create the brain mask with grey matter and white matter

clear Vi Vo
Vi = strvcat([subjdir 'c1T1F.nii'], [subjdir 'c2T1F.nii']); % grey matter (c1T1F.nii) and white matter (c2T1F.nii) tissue maps for brain mask
Vo = [subjdir 'GMWM_T1F_mask.nii']; % The brain mask 
spm_imcalc(Vi, Vo, '(i1+i2)>0.9'); % Creating the brain mask by adding grey matter and white matter

% Dilate and erode the brain mask using spherical neighborhood to fill in
% all the voids in the brain mask which may be due to tissue
% misclassification from segmentation
V=spm_vol(Vo);
c=spm_read_vols(V); % Load the brain mask
[xx,yy,zz] = ndgrid(-5:5);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 4.0; % define a three-dimensional spherical neighborhood of 2 unit radius 
c=imdilate(c,nhood);
c=imdilate(c,nhood);
c = imerode(c,nhood);
c = imerode(c,nhood);
V.fname='GMWM_T1F_mask.nii';
spm_write_vol(V,c); % save the smoothed brain mask

%%  Coregister the T1 MRI (along with tissue probability maps and brain masks) to the post-implantation CT. 

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[subjdir 'CTpost.nii,1']}; % post implantation CT
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[subjdir 'mT1F.nii,1']}; % T1 MRI
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[subjdir 'c1T1F.nii,1'] % grey matter
                                                   [subjdir 'c2T1F.nii,1'] % white matter
                                                   [subjdir 'c3T1F.nii,1'] % CSF
                                                   [subjdir 'brain_T1F_mask.nii,1'] % unsmoothed brain mask with CSF
                                                   [subjdir 'brain_T1F_mask_ero.nii,1'] % smoothed brain mask with CSF
                                                   [subjdir 'GMWM_T1F_mask.nii,1'] % smoothed brain mask without CSF
                                                   };
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r'; % save all the tissue probability maps and brain masks with prefix 'r'
matlabbatch{2}.spm.spatial.coreg.write.ref = {[subjdir 'CTpost.nii,1']};
matlabbatch{2}.spm.spatial.coreg.write.source = {[subjdir 'mT1F.nii,1']};
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r'; % save coregistered T1 MRI as 'rmT1F.nii'
spm_jobman('run', matlabbatch);
clear matlabbatch

%% Create brain masked T1 MRI
Vi = strvcat([subjdir 'rmT1F.nii'], [subjdir 'rbrain_T1F_mask_ero.nii']); % use smoothed brain mask with CSF
Vo = [subjdir 'brain_rmT1F.nii']; % Define the brain mased T1 MRI
spm_imcalc(Vi, Vo, 'i1.*i2'); % save brain masked T1 MRI
clear Vi Vo V

    %% Create CT difference image 
    filerCT1 = [subjdir '/rCTpre.nii']; % Load pre-impantation CT
    fileCT2 = [subjdir 'CTpost.nii']; % Load post-impantation CT
    VolCT2 = spm_vol(fileCT2); % reference image
    Vi = strvcat(fileCT2, filerCT1);
    Vo = [subjdir 'CTdiff.nii']; % The CT difference image
    spm_imcalc(Vi, Vo, 'i1-i2'); % Subtracting the pre-implantation CT from post-implantation CT
    
    %% Apply the brain mask to CT difference image
    clear Vi Vo
    Vi = strvcat([subjdir 'CTdiff.nii'], [subjdir 'rbrain_T1F_mask_ero.nii']); % use smoothed brain mask with CSF
    Vo = [subjdir 'brain_masked_CTdiff.nii']; % Brain masked CT difference image
    spm_imcalc(Vi, Vo, 'i1.*(i2)'); % Save brain masked CT difference image
