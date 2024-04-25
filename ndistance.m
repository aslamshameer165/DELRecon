function dist=ndistance(crd,clusters,weight)

% traditional K-means distance
mode='normal'; % normal mode for unweighted Euclidean distance 
% Choose weighted mode for weighted distance

switch mode
    case 'normal'
        
        for i=1:length(crd)
            dist(i)=norm(crd(i,:)-crd(clusters(i),:)); % Distance of all data points from current cluster
        end
        
    case 'weighted'
        for i=1:length(crd)
            dist(i)=weight(i).^p * norm(crd(i,:)-crd(clusters(i),:),p); % Weighted distance of all data points from current cluster
        end
end