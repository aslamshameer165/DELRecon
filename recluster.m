function [dist,clusters] = recluster(crd,dist,clusters,weight,i)

%RECLUSTER  [dist,clusters] = recluster(crd,dist,clusters,weight,i)
%       Reassign data points to the hubs
% i index to new cluster hub

n = size(crd,1);
temp = i*ones(n,1);
newdist=ndistance(crd,temp,weight);  % Calculate the distance of all points to new hub

    for r=1:n
        if newdist(r) < dist(r)     %Reassign the cluster if the new distance is less
           dist(r)=newdist(r); % if new distance is less than the distance to previous cluster, assign it to new cluster
           clusters(r)=i; %otherwise leave it with previous cluster
        end     

     end        

end
