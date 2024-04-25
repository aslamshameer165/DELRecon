function resDist=svdfitDist(a,n,P)
% a: center of cluster
% n: direction vector of SVD line
% P: Data ponts

for i=1:size(P,1)
   resDist(i,1)=norm((a-P(i,:))-((a-P(i,:))*n)*n') ; % Euclidean distance of P from the cluster(gives column vector)
end