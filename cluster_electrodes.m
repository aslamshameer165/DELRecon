%% Cluster Electrodes
% This script implement the Module 2 (Electrode Detection) and Module 3
% (Individual Contact Labeling) for DELRecon
%Performs clustering of sEEG electrodes. Thresholds CT image to
%leave only voxels corresponding to electrodes and clusters outermost
%voxels to identify individual handles. Plots interactive figure for visual
%inspection of initial clustering accuracy. If misclassification occurs at
%this stage, you can use the brush tool to draw boxes around electrode
%handles to specify CLUSTER or NOCLUSTER.
%
%Requirements/dependencies: Smoothed brain mask with CSF from T1 MRI and CT
%difference image (both can be obtained using the function CT_MR_preprocess
%if they do not already exist - naming convention must be as follows:
%'brain_masked_CTdiff.nii' 'brain_rT1F_mask_ero.nii').
%
%Inputs (line 31-40): sub, nEle, firstThreshold
%where sub is a char array containing patient ID, nEle is the number of
%electrodes for the patient and firstThreshold is the initial thresholding value
%for the CT image
%
%Output of initial clustering is saved by default to 'init_CL.mat', and
%contact localisation to 'Contacts.mat', both in working directory (i.e. 'basedir')


clear all
close all
clc

%%  Define path to the data and dependent scripts
%addpath(genpath('E:\PROJECT\SEEG segmentation\Edited\scripts\iElectrodes')); % Path to the iElectrodes script folder
addpath(genpath('E:\PROJECT\SEEG segmentation\Edited\DELRecon')); % Path to the main script folder

%Inputs - change per patient (Keep default treshold value in most cases unless required to change)
sub = 'subject'; % subject folder to be processed
nEle = 7; % Define number of electrodes implanted
firstThreshold = 650; % Define the threshold to create the electrode image from CT difference image
clr = lines(nEle); % for color of different lines based on number of electrodes
basedir = strcat('E:\PROJECT\SEEG segmentation\Edited\DELRecon\',sub,'\Imaging\'); %path to the data to be processed
cd(basedir);

%% Load brain masked CT difference image
fprintf('Loading CT data \n');
CT_masked = fullfile(basedir,'brain_masked_CTdiff.nii'); % the brain masked CT difference image
V=spm_vol(CT_masked);
ima=spm_read_vols(V);
%% Load brain mask
fprintf('Loading brain mask \n');
brain_mask = fullfile(basedir,'rbrain_T1F_mask_ero.nii'); % the binary brain mask
V_mask=spm_vol(brain_mask);
ima_mask=spm_read_vols(V_mask);

%% Erode the brain mask to get the seed mask close to cortical surface
fprintf('Eroding brain mask \n');
[xx,yy,zz] = ndgrid(-5:5);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 4.0; % define a three-dimensional spherical neighborhood of 2 unit radius
ima_mask_erode2 = imerode(ima_mask,nhood); % First stage erossion
ima_mask_erode2 = imerode(ima_mask_erode2,nhood); % Second stage erossion
ima_mask_erode=ima_mask_erode2;
ima_mask_seed=ima_mask-ima_mask_erode2; % Subtract the erodded brain mask from the original brain mask to generate a small strip of cortical surface mask

%% binarize CT data with a given threshold
fprintf('Binarising CT using given threshold \n');
ima_bw=(ima>firstThreshold); % Value 1 for the voxels of CT difference image higher than the threshold (Electrode image mask)
ind_seed = find (ima_mask_seed & ima_bw); % Electrode voxels in the cortical surface strip mask (Electrode handles)
[x_seed, y_seed, z_seed] = ind2sub(size(ima), ind_seed); % Coordinates of electrode handles within the cortical surface strip mask
[x, y, z] = ind2sub(size(ima), find(ima_mask>0)); % Coordinates of brain mask
[x_all, y_all, z_all] = ind2sub(size(ima), find(ima_bw & ima_mask_erode)); % coordinates of electrode voxels beyond the cortical surface strip mask

% Keep only tip of electode handles
validIndx = find(x_seed<(min(x)+(max(x)-min(x))*0.25) | x_seed>(max(x)-(max(x)-min(x))*0.25) | z_seed>(min(z)+(max(z)-min(z))*0.5)); % Voxels of electrode handles close to the cortical surface

% save some memory
clear ima_mask_erode ima_mask_seed;

%% Initial clustering to determine the commencement of each electrode at cortical surface
if ~isempty(dir(fullfile(basedir,'initCL.mat'))) % check if the initial clustering has already being performed
    load(fullfile(basedir,'initCL.mat')); % load if initCL.mat(initial clustering result) exist
else
   % check initial clusters, correct if problem occurs here using GUI 1
   
   % Cluster the tips of electrode handle into number of electrodes.
   % Use the function for customized clustering from iElectrodes toolbox to give fixed solution.
   % clustering(A,B,C): A) Number of clusters (Electrodes), B) Coordinates of data points (Tips of electrode handles), C) Weight (Unweighted clustering) 
   % Clustering results in GS_seed: center of mass of each cluster and clusters_seed: cluster label corresponding to each data point
   [GS_seed,clusters_seed]=  clustering (nEle,[x_seed(validIndx) y_seed(validIndx) z_seed(validIndx)],ones(length(x_seed(validIndx)),1));
    
   % Current Cluster
    x_clustered=x_seed(validIndx); % x-coordinate of tip of electrode handle
    y_clustered=y_seed(validIndx); % y-coordinate of tip of electrode handle
    z_clustered=z_seed(validIndx); % z-coordinate of tip of electrode handle
    clustered=clusters_seed; % cluster label corresponding to each data point
    
    % Quality check using GUI 1. Do we have all electrodes in separate clusteres?
    f=figure(1); % open figure
    f.WindowState='maximized'; %maximize the figure window
    hold on;
    clr = lines(nEle); % for color of different lines
    set(f,'Position',[400 80 1100 700]); % width and height of figure
    S.seed=scatter3(x_clustered, y_clustered, z_clustered,36,clr(clustered,:),'Marker','.'); % Plot tip of electrode handles
    S.GS=scatter3(GS_seed(:,1),GS_seed(:,2),GS_seed(:,3),100,clr(1:nEle,:),'Marker','o','LineWidth',3); % Plot center of mass of each cluster
    S.all=scatter3(x_all, y_all, z_all,1,'Marker','o'); %plots remaining electrode voxels beyond the cortical surface strip 
    view(3);axis vis3d, box on;rotate3d on;
    xlabel('x'),ylabel('y'),zlabel('z');
    title('If needed, get CLUSTER and NOCLUSTER variables using the brush tool');
    
    % Referesh button to update the clustering results with the user defined clusters.
    p = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.85 .05 .1 .08],...
        'String','Refresh',...
        'CallBack','uptFigure');
    
    % Plot or remove the electrode voxels beyond the cortical surface strip from the figure
    p2 = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.05 .85 .1 .08],...
        'String','Vox on/off',...
        'CallBack','uptFigure2');
    
    % Plot or remove the tip of electrode handles from the figure
    p3 = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.05 .75 .1 .08],...
        'String','Seed on/off',...
        'CallBack','uptFigure3');
    
    % Plot or remove the center of mass of each cluster from the figure
    p4 = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.05 .65 .1 .08],...
        'String','GS on/off',...
        'CallBack','uptFigure4');
    
    % To close the figure
    q = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.05 .05 .1 .08],...
        'String','Close',...
        'CallBack','close(gcbf)');
    
    % To discard the user defined variables for clustering
    q2 = uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.85 .85 .1 .08],...
        'String','Restore',...
        'CallBack','clearCLUSTERs');
    waitfor(f); % wait until figure is closed
    
    % save the initial clustering results as initCL.mat in the working directory
    save(fullfile(basedir,'initCL.mat'),'clusters_seed','x_clustered','y_clustered','z_clustered','GS_seed','clustered');
    
end




%% For each cluster, define a SVD line
ima_clustered=zeros(size(ima));  % Define a blank image
% set clusters from seeds
for i=1:nEle %to save clustered voxels of electrodes in a image file as the label of clusters
    ima_clustered(sub2ind(size(ima_clustered),x_clustered(clustered==i),y_clustered(clustered==i),z_clustered(clustered==i))) = i;
end
% the rest voxels to be clustered within Erode2 Masks
ind_rest = find (ima_mask_erode2 & ima_bw); % electrode in the region beyond the cortical surface strip
[x, y, z] = ind2sub(size(ima), ind_rest); % coordinates of electrode voxels beyond the cortical surface strip
newCor = setdiff([x y z],[x_clustered y_clustered z_clustered],'rows'); %electrode voxels which are not clustered
%coordinates of unclustered electrodes beyond the cortical surface strip
x=newCor(:,1);
y=newCor(:,2);
z=newCor(:,3);

% Fit SVD lines per electrode
for i=1:nEle
    % fit a line
    D=[x_clustered(clustered==i) y_clustered(clustered==i) z_clustered(clustered==i)];% Coordinated of ith cluster
    % svdfit implements the Singular Value Decomposition to fit a line on the data points
    [p0(i,:),d(:,i)]=svdfit(D);% Define the svd line for ith cluster having coordinates in D
    % p0=center of the cluster (coordinates in ith row for cluster i)
    % d=direction vector of the cluster (vector in ith column for cluster i)
    % average is arranged in single column
    % direction is arranged vertically in each column
    maxDist(i) = 20;  % Define maxDist of each cluster with the value 20
end

%% Stepwise clustering of remaining electrode voxels
allSVDdist=Inf; % initialize the maximum distance of any point to its closest cluster
while (1)
    clear allDist;
    for i=1:nEle % Calculate the  Euclidean distance of remaining electrode voxels from each cluster
        % svdfitDist calculates the  Euclidean  distance of the data point from the p0 along d
        allDist(:,i)=svdfitDist(p0(i,:),d(:,i),[x_all y_all z_all]); %ith column has distance of remaining points from ith cluster
    end
    [~, tempCluster]=min(allDist,[],2); %tempCluster= cluster closest to each point
    allSVDdist_cur=max(allDist(sub2ind(size(allDist),1:size(allDist,1),tempCluster'))); % the maximum distance of any point to its closest cluster
    
    % stop tracking if the SVD distance is stable(i.e., the maximum distance of any point to its closest cluster reduced by not more than 2 units)
    if (allSVDdist - allSVDdist_cur ) < 2
        break;
    % Otherwise, continue tracking
    else
        allSVDdist=allSVDdist_cur; % Update the new maximum distance of any point to its closest cluster
        disp(['Current max of SVD distance: ' num2str(allSVDdist)]);
        clear allDist;
        steps = 100; % Define the number of data points to be tracked at a time
        k=0; 
        % Run the following loop 20 times or until any data point remains unclustered
        while ~isempty(x) && k<20
            disp([num2str(length(x)) ' voxels remaining ...']); % Number of unclustered voxels
            stepClusters=min(length(x),steps); % number of data points to be clustered (maximum 100 or whichever remained)
            clear allDist allIndx allDistfromGS;
            for i=1:nEle
                % Distance of each data point from the center of mass of ith cluster
                allDistfromGS(:,i)=sqrt(sum(([x y z]-repmat(GS_seed(i,:),[size(x),1])).^2,2));
                % Redefine SVD lines for clustered data points
                D=[x_clustered(clustered==i) y_clustered(clustered==i) z_clustered(clustered==i)];
                [p0(i,:),d(:,i)]=svdfit(D);
                
                allDist(:,i)=svdfitDist(p0(i,:),d(:,i),[x y z]); % distance of each data point from ith SVD line
            end
            
            % Find voxels closest to the center of mass of the clusters
            CostFunc=allDistfromGS; % Define cost function
            [tempVal,tempIndx]=sort(CostFunc(:)); % arrange the distances in ascending order
            
            % Choose maximum 100 data points or whichever remained
            % unclustered based on the distance from the center of mass of clusters
            [tempVox, tempCluster]=ind2sub(size(CostFunc), tempIndx(1:stepClusters));
            
            % Define distance of the choosen data points from its closest
            % cluster along the corresponding SVD line
            tempDist = allDist(sub2ind(size(allDist),tempVox,tempCluster));
            
            strangIndx=tempDist>maxDist(tempCluster)'; % To discard the data points if its distance from closest cluster is more than the maxDist(i.e., 20)
            tempVox(strangIndx)=[];
            tempCluster(strangIndx)=[];
            
            if length(unique(tempVox))~=stepClusters % Warn in case of overlapped voxels being clustered
                warning('overlapping clusters');
            end
            
            for i=1:length(tempVox)%Update the cluster image based on the new data points being clustered
                ima_clustered(sub2ind(size(ima_clustered),x(tempVox(i)),y(tempVox(i)),z(tempVox(i)))) = tempCluster(i);
            end
            
            % Add new clusteres to already defined clusters
            clustered = [clustered; tempCluster];
            x_clustered = [x_clustered;x(tempVox)];
            y_clustered = [y_clustered;y(tempVox)];
            z_clustered = [z_clustered;z(tempVox)];
            

            % Remove the clustered data points from the list of unclustered
            x(tempVox)=[];
            y(tempVox)=[];
            z(tempVox)=[];
            
            % Remove the discarded data points (when their distance from closest cluster was more than the maxDist)
            x(strangIndx)=[];
            y(strangIndx)=[];
            z(strangIndx)=[];
            k=k+1; % Run for 20 times
        end
    end
end

%% Further tracking for the rest of the voxels
% fit the remaining voxels along the SVD line
clear allDist;
for i=1:nEle % To calculate the distance of unclustered data points from the p0 along d
    allDist(:,i)=svdfitDist(p0(i,:),d(:,i),[x y z]);
end
[~, tempCluster]=min(allDist,[],2); % Determine the cluster closest to each point, hence, define their cluster label
for i=1:length(x) %Update the cluster image based on the new data points being clustered 
    ima_clustered(sub2ind(size(ima_clustered),x,y,z)) = tempCluster; 
end

% Add new clusteres to already defined clusters
clustered = [clustered; tempCluster];
x_clustered=[x_clustered;x];
y_clustered=[y_clustered;y];
z_clustered=[z_clustered;z];



% Redefine the SVD line
for i=1:nEle
    D=[x_clustered(clustered==i) y_clustered(clustered==i) z_clustered(clustered==i)];
    [p0(i,:),d(:,i)]=svdfit(D);
    
end



  

%% Recluster based on the final SVD line
clear allDist;
clear allDistfromMean;
clear allDistfromGS;
clear x y z;

[x, y, z] = ind2sub(size(ima), find(ima_bw));%% Coordinates of all electrode voxels

for i=1:nEle
    % The distance of data points from the center of SVD line along the line
    allDist(:,i)=svdfitDist(p0(i,:),d(:,i),[x y z]);
    % Euclidean distance to the center of the line
    allDistfromMean(:,i)=sqrt(sum((repmat(p0(i,:),length(x),1)-[x y z]).^2,2));
    % Euclidean distance to the center of mass of initial clusters
    allDistfromGS(:,i)=sqrt(sum(([x y z]-repmat(GS_seed(i,:),[size(x),1])).^2,2));
end

% Distance metric to define final clusters
% = combined euclidean from the SVD line along their direction, center of SVD line and the initial clusters
[~, tempCluster]=min(allDist.*20+allDistfromMean+allDistfromGS,[],2);

% Update all the clusters based on the distance metric
clustered=tempCluster;
x_clustered=x;
y_clustered=y;
z_clustered=z;

% Refit the final SVD line
for i=1:nEle
    D=[x_clustered(tempCluster==i) y_clustered(tempCluster==i) z_clustered(tempCluster==i)];
    [p0(i,:),d(:,i)]=svdfit(D);
end

    

% Update the final cluster image
ima_clustered=zeros(size(ima));
for i=1:length(x_clustered)
    ima_clustered(sub2ind(size(ima_clustered),x_clustered,y_clustered,z_clustered)) = tempCluster;
end

%% Plot brain surface and electode handles
% Recheck if all electrodes are clustered correctly
% Based on this plot, we can assign electrode labels to clusteres
% Use the plot to determine the clusterID of each electrode
% Put result in a text file "EleMap" in a format (%ElectodeName %contactNumber %clusterID)
% The electrode description in EleMap should be arranged in the order of their clusterID
% Example:
% P 14 1
% Y' 16 2
% Y 16 3
% P' 12 4
% X 16 5
% O' 12 6
% ...

% To plot the 3D brain with all data points and SVD lines
h=plotBrainandElec(ima_mask_erode2,ima_clustered,x,y,z,[],[],[],clustered,GS_seed,clr,1);
t=-100:100; % Define the length of SVD line to be plotted
for k=1:nEle
    for i=1:length(t)
        P(i,:) = p0(k,:) + d(:,k)'.*t(i); % Determine the SVD line to be plotted
    end
    % Plot the data points with the defined color
    scatter3(h.ax2,P(:,1),P(:,2),P(:,3),5,repmat(clr(k,:),size(P,1),1),'.');
    scatter3(h.ax1,P(:,1),P(:,2),P(:,3),5,repmat(clr(k,:),size(P,1),1),'.');
    clear P;
end




%% Individual Contact Labeling
% Read contact file in a format (%ElectodeName %contactNumber %clusterID)
% Can also pass the 4th and 5th columns, which limit the range [4th column, 5th column] of the SVD
% Can also pass the 6th and 7th columns, which limit the range [6th column, 7th column] of the voxel intensity
% line, to use for exclusively masking heads and tails.

clear ima_bw ima_mask_erode2 ima_mask_erode3
fid = fopen(fullfile(basedir,'EleMap.txt')); %Open a text file
FC = textscan(fid, '%s%f%f%f%f%f%f'); % Scan the content of the file along the row
fclose(fid); 
Contacts.contName=FC{1}; % 1st column of file defines the name of electrode
Contacts.contNum=FC{2}; % 2nd column defines the number of contacts in each electrode
Contacts.cluster=FC{3}; % 3rd column defines the cluster ID of each electrode
Contacts.xrange = [FC{4} FC{5}]; % 4th and 5th column defines the length of the electrode
Contacts.yrange = [FC{6} FC{7}]; % 6th and 7th column defines the intensity range
Contacts.p0=p0; % Center of SVD lines
Contacts.d=d; % Direction vector of SVD lines

% Iterate through each electrode
for i=1:nEle
pEle=i;
    clear f S allP P bbb linesp linedist indx xx GS2 clusters2 p1 p2 p3 p4 q1 q2 q3 q4
    f=figure(pEle); % Open a new figure for each electrode
    f.WindowState='maximized';
    hold on;
    set(f,'Position',[400 80 1100 700]);
    allP=[x_clustered(clustered==pEle),y_clustered(clustered==pEle),z_clustered(clustered==pEle)]; % Coordinates of data points in ith electrode
    % only consider voxels with intensities more than firstThreshold
    P=allP(ima(sub2ind(size(ima),allP(:,1),allP(:,2),allP(:,3)))>firstThreshold,:); % Coordinates after thresholding
    
    bbb=ima(sub2ind(size(ima),P(:,1),P(:,2),P(:,3))); % Image of ith electrode after thresholding
    linesp=svd3D_1D(p0(pEle,:), d(:,pEle)',P);% Calculate the unit vector from SVD center point to the P along the line (This is to projest the electrodes to 1D)
    linedist=svdfitDist(p0(pEle,:),d(:,pEle),P);% Distance of P from SVD center point along the line
    
    uu=0; % Initialize scaling parameter for figure zoom-in/out
    
    % Plot the data points in 1D along the SVD line
    S.sc=scatter(linesp(linedist<20),bbb(linedist<20)-uu,1,repmat(linedist(linedist<20)./max(linedist),1,3));
    ylim([uu,max(bbb(linedist<20))-uu])
    if isnan(Contacts.xrange(pEle,1)) % Define the lower range of electrode length
        Contacts.xrange(pEle,1)=min(linesp);
    end
    if isnan(Contacts.xrange(pEle,2)) % Define upper range of electrode length
        Contacts.xrange(pEle,2)=max(linesp);
    end
    if isnan(Contacts.yrange(pEle,1)) % Define lower range of intensity
        Contacts.yrange(pEle,1)=min(bbb(linedist<20));
    end
    if isnan(Contacts.yrange(pEle,2)) % Define upper range of intensity
        Contacts.yrange(pEle,2)=max(bbb(linedist<20));
    end
    
    % clustering to get individual contact
    % Define the data points with in the defined ranges and distance of less then 2 units from SVD line
    indx=linedist<2 & linesp<=Contacts.xrange(pEle,2) & linesp>=Contacts.xrange(pEle,1) & bbb(linedist<20)<=Contacts.yrange(pEle,2) & bbb(linedist<20)>=Contacts.yrange(pEle,1);
    
    % Coordinates of the data points under above constraints
    xx=sub2ind(size(ima),P(indx,1),P(indx,2),P(indx,3));
    
    % Determine the clusteres in the data points under constraints to define the contacts
    [GS2,clusters2]=  clustering(Contacts.contNum(pEle),linesp(indx),ima(xx));
    
    linespGS = svd3D_1D(p0(pEle,:), d(:,pEle)',GS_seed(pEle,:));
%   S.pt=scatter(GS2(:,1),zeros(Contacts.contNum(pEle),1)+max(bbb(linedist<20))-uu,'r');
    % Plot the contacts
    S.pt=scatter(GS2(:,1),zeros(Contacts.contNum(pEle),1)+(Contacts.yrange(pEle,1)+Contacts.yrange(pEle,2))/2,'r');
    % Plot the lower range of electrode length
    S.ln1=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,1)],[0 max(bbb(linedist<20))],'r');
    % Plot the upper range of electrode length
    S.ln2=plot([Contacts.xrange(pEle,2) Contacts.xrange(pEle,2)],[0 max(bbb(linedist<20))],'r');
    % Plot lower intensity range
    S.ln1y=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,2)],[Contacts.yrange(pEle,1) Contacts.yrange(pEle,1)],'b');
    % Plot upper intensity range
    S.ln2y=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,2)],[Contacts.yrange(pEle,2) Contacts.yrange(pEle,2)],'b');


% Displays the number of contacts    
static_text1 = uicontrol('Style', 'text',...
    'Units', 'Normalized',...
    'Position', [0.52, 0.95, 0.13, 0.03],...
    'String', 'No. of contacts', 'FontSize', 12, ... % Set the font size
    'FontWeight', 'bold');
text_box1 = uicontrol('Style', 'edit',...
    'Units', 'Normalized',...
    'Position', [0.7, 0.95, 0.03, 0.03],...
    'String', num2str(Contacts.contNum(pEle)),...
    'Tag', 'Contacts',... % Set the Tag property
    'Enable', 'off',...
    'Callback', 'disp(get(gcbo, ''String''))');

% Displays electrode number
static_text2 = uicontrol('Style', 'text',...
    'Units', 'Normalized',...
    'Position', [0.22, 0.95, 0.08, 0.03],...
    'String', 'Electrode', 'FontSize', 12, ... % Set the font size
    'FontWeight', 'bold');
text_box2 = uicontrol('Style', 'edit',...
    'Units', 'Normalized',...
    'Position', [0.32, 0.95, 0.03, 0.03],...
    'String', num2str(i),...
    'Tag', 'Electrode',... % Set the Tag property
    'Enable', 'off',...
    'Callback', 'disp(get(gcbo, ''String''))');

% Define pannel for zoom button
zoomPanel = uipanel('Title', 'Zoom', 'Position', [0.01, 0.8, 0.08, 0.15]);

% Add buttons to the 'zoom' panel
% To zoom-in the plot
q3 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Parent', zoomPanel,...
    'String', 'IN',...
    'Position', [0.1, 0.65, 0.8, 0.3],...
    'CallBack', 'inup');
% To zoom out the plot
q4 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Parent', zoomPanel,...
    'String', 'OUT',...
    'Position', [0.1, 0.35, 0.8, 0.3],...
    'CallBack', 'indw');

% Define panel for intensity ranges
panelIntensityRange = uipanel('Title', 'Intensity range', 'Position', [0.91, 0.6, 0.08, 0.15]);
% To change lower intensity range
py3 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Position', [0.05, 0.7, 0.9, 0.25],... % Adjust the position within the panel
    'String', 'Lower',...
    'CallBack', 'yrange1',...
    'Parent', panelIntensityRange);
% To change upper intensity range
py4 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Position', [0.05, 0.4, 0.9, 0.25],... % Adjust the position within the panel
    'String', 'Upper',...
    'CallBack', 'yrange2',...
    'Parent', panelIntensityRange);

% Panel to manipulate the number of contacts
changeContact = uipanel('Title', 'Contacts', 'Position', [0.91, 0.4, 0.08, 0.15]);
% To add new contact
p1 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Position', [0.05, 0.7, 0.9, 0.25],... % Adjust the position within the panel
    'String', 'Add',...
    'CallBack', 'adcontact',...
    'Parent', changeContact);
% To remove a contact
p2 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Position', [0.05, 0.4, 0.9, 0.25],... % Adjust the position within the panel
    'String', 'Remove',...
    'CallBack', 'rmcontact',...
    'Parent', changeContact);

% Panel to redefine electrode length
panelElectrodeLengthRange = uipanel('Title', 'Electrode length range', 'Position', [0.91, 0.8, 0.08, 0.15]);
% To change lower range of electrode length
p3 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Position', [0.05, 0.7, 0.9, 0.25],... % Adjust the position within the panel
    'String', 'Lower',...
    'CallBack', 'range1',...
    'Parent', panelElectrodeLengthRange);
% To change upper range of electrode length
p4 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Position', [0.05, 0.4, 0.9, 0.25],... % Adjust the position within the panel
    'String', 'Upper',...
    'CallBack', 'range2',...
    'Parent', panelElectrodeLengthRange);

% To close the figure
q1 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Position', [0.01, 0.05, 0.08, 0.08],...
    'String', 'Close',...
    'CallBack', 'close(gcf)');

% To reset to the default parameters
q2 = uicontrol('Style', 'PushButton',...
    'Units', 'Normalized',...
    'Position', [0.01, 0.15, 0.08, 0.08],...
    'String', 'Reset',...
    'CallBack', 'rst');

    waitfor(f);
    

    % coordinate for each contact, the first index always the one
    % furthest to the handle
    % Sort the contacts based on the direction of implantation
    GScontacts = sort(GS2(:,1),'ascend');
    if dist(linespGS,GScontacts(1))<dist(linespGS,GScontacts(end))
        GScontacts = sort(GS2(:,1),'descend');
    end
    
    % Update the EleMap.txt file based on the determined electrode contacts
    for j=1:length(GScontacts)
        Contacts.contcoor{pEle}(j).GS=GScontacts(j);
        tempcoor=round(svd1D_3d(p0(pEle,:), d(:,pEle)',GScontacts(j)));
        % find the max intensity and its coordinates in a 3x3x3 cube
        bnd = 1;
        [x y z] = meshgrid((tempcoor(1)-bnd):(tempcoor(1)+bnd),(tempcoor(2)-bnd):(tempcoor(2)+bnd),(tempcoor(3)-bnd):(tempcoor(3)+bnd));
        [~, b]=max(ima(sub2ind(size(ima),x(:), y(:) ,z(:))));
        Contacts.contcoor{pEle}(j).coor=[x(b), y(b) ,z(b)];
    end
     Contacts.contNum(i)=length(GS2);
end

% Update the EleMap.txt file based on the user defined parameters
fid = fopen(fullfile(basedir,'EleMap.txt'),'r');
FC = textscan(fid, '%s%f%f%f%f');
fclose(fid);
fid = fopen('EleMap_NEW_updated.txt', 'w');
for i=1:nEle
FC{2}(i)=Contacts.contNum(i);  
FC{4}(i)=Contacts.xrange(i,1);
FC{5}(i)=Contacts.xrange(i,2);

fprintf(fid, '%s\t%f\t%f\t%f\t%f\n', FC{1}{i}, FC{2}(i), FC{3}(i), FC{4}(i), FC{5}(i));
end
fclose(fid);

save(fullfile(basedir,'Contacts'),'Contacts');




%%   Create electrode image
CT_masked = fullfile(basedir,'rCTpre.nii'); % Read a sample image
V=spm_vol(CT_masked);
ima_mask=spm_read_vols(V);
% s=uint16(zeros(size(ima)));
s=double(zeros(size(ima_mask))); % Define a blank image based on the resolution of sample image
clear k h

[xx,yy,zz] = ndgrid(-2:2);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 2; % create the spheres at the position of each contact. Change the radius of the displayed contacts here.

% Define spheres at the position of each contact to create a electrode image
for i=1:length(Contacts.contNum)
    for j=1:length(Contacts.contcoor{1,i})
        k=Contacts.contcoor{1,i}(j).coor;
        % The value at each sphere defines the position of the contact in the image
        % The ones and tens place digit defines the contact number
        % The digits preceeding the contacts digits defines the electrode number
        s(k(1)-2:k(1)+2,k(2)-2:k(2)+2,k(3)-2:k(3)+2)=(nhood*i*100)+(nhood*j); 
    end
end
    
nnn=V;

nnn.fname=[pwd '\all_electrodes.nii'];
spm_write_vol(nnn,s); % Save the electrode image

