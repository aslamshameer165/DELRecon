%% To create an ITKSnap software compatible label file for electrode image

clear
clc
%% Define
sub = 'subject'; % subject folder to be processed
nEle = 7; % Define number of electrodes implanted
basedir = strcat('E:\PROJECT\Swansea\SEEG segmentation\Edited - OneDrive_1_10-6-2023\DELRecon\',sub,'\Imaging\'); %path to the data to be processed
cd(basedir);

%% Create a label file
fid = fopen('Contact_labels.txt', 'w'); % open a file to define contact labels

% Define all the columns compatible for ITKSnap software
FC=cell(1,8);
load('Contacts.mat'); % Load the contact information
a1=[];
a2=[];
a3=[];
a4=[];
a5=[];
a6=[];
a7=[];
a8=[];
% FC = textscan(fid, '%d%d%d%d%d%d%d%s');

% Load the tissue probability maps
GM = [pwd '\rc1T1F.nii']; % Grey matter
WM = [pwd '\rc2T1F.nii']; % White matter
CSF = [pwd '\rc3T1F.nii']; % CSF
V=spm_vol(GM);
GM=spm_vol(GM);
WM=spm_vol(WM);
CSF=spm_vol(CSF);
GM=spm_read_vols(GM);
WM=spm_read_vols(WM);
CSF=spm_read_vols(CSF);

% Write the label file
for i=1:nEle
    for j=1:Contacts.contNum(i)
        s=Contacts.contName{i};
        a1=[a1 (i*100)+j]; % Contact description
        % Fill in the constants to define the labels
        a2=[a2 255];
        a3=[a3 0];
        a4=[a4 0];
        a5=[a5 1];
        a6=[a6 1];
        a7=[a7 1];
        
        c=Contacts.contcoor{1,i}(j).coor; % Get the centroid of the contact
        % Classify the contact into one of the tissue types
        if GM(c(1),c(2),c(3))>0.5
            s=[s '-GM'];
        elseif WM(c(1),c(2),c(3))>0.5
            s=[s '-WM'];
        elseif CSF(c(1),c(2),c(3))>0.5
            s=[s '-CSF'];
        else
            s=[s '-NA'];
        end
        s=['"' s '"'];
        a8=[a8 string(s)];
    end
end

% Define a clear label
t=[];
t=[t "0"];
t=[t "0"];
t=[t "0"];
t=[t "0"];
t=[t "0"];
t=[t "0"];
t=[t "0"];
t=[t,"Clear Label"];
t(8)=strcat('"',t(8),'"');
a=[num2str(a1') num2str(a2') num2str(a3') num2str(a4') num2str(a5') num2str(a6') num2str(a7') a8'];
a=[t;a];
fid = fopen('Contact_labels.txt', 'w');
for i=1:length(a) % write the label descriptions to the file
fprintf(fid, '\t\t%s\t\t\t%s\t\t%s\t\t%s\t\t\t\t\t\t\t\t%s\t\t%s\t\t%s\t\t\t\t%s\n', a(i,1), a(i,2), a(i,3), a(i,4), a(i,5), a(i,6), a(i,7), a(i,8));
end
fclose(fid);