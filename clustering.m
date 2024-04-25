function [GS,clusters]=clustering (nClusters,crd,weight)
% This is a cutomized clustering algorithm to provide the consistant clustering results

% function GS=clustering (nClusters,crd)
% compute nClusters using k-means algorithm. 

% function GS=clustering (nClusters,crd,weight)
% calculate cluster coordinates using a weigthed center of mass of each cluster

% GS: center of mass of each cluster
% clusters: cluster labels corresponding to each data point

l=size(crd,1); % define an array of size of number of data points

% Normalize the weights
temp=weight-min(weight);
normWeight=temp/max(weight);


% Initialize all the data points to cluster 1
clusters=ones(l,1);

hubs=zeros(nClusters,1);    % total number of clusters, initially 0 in all
hubs(1,1)=1;                % start with 1 cluster

dist = ndistance(crd,clusters,normWeight); %distance of all crd to crd(1)
                 %hubs
counter =1;      %number of hubs
continuar = 1;   % 1=true 0=false indicates whether to continue forming new cluster
                 


while continuar
    counter = counter + 1;  % adding new hub
    [m,i]=max(dist);        % m = Distance of farthest data point from the last cluster
                            % i = Index of farthest data point from last cluster
      
    hubs(counter)=i;        % assign farthest data point from last cluster as another hub (another cluster)
 
    [dist,clusters] = recluster(crd,dist,clusters,normWeight,i);  % Clusters points to nearest hub and recalculate the distance
                                                                  
                                                            
%    maxdist=max(dist);      %returns distance of point farthest from its hub
%    continuar = farout(counter,hubs,crd,maxdist); %Checks the stop condition
% si la distancia es mayor a 8mm continuo
%     if maxdist>10

% use known number of cluster to stop
     if counter<nClusters % Stop clustering if all data points are exhausted
        continuar=1;
     else
         continuar=0;
     end
    
end

%% calculate final locations for each cluster

GS=zeros(counter,3); % Initialize a variable to hold center of mass of each cluster

% Calculate the average coordinate of each cluster
% weighted center of mass
if nargin==3 
    for i=1:counter
        weight=double(weight);
        ind=clusters==hubs(i);
        W=sum(weight(ind));
        GS(i,:)=sum(diag(weight(ind))* crd(ind,:),1)*(1/W);
    end

else
% non-weigthed center of mass    
    for i=1:counter
        GS(i,:)=mean(crd(clusters==hubs(i),:),1);
    end
end

%% reorganize output - Number clusters 1:N 
if nargout>1
    temp=clusters;
    for i=1:counter
        clusters(temp==hubs(i))=i; % Define cluster labels
    end
end

