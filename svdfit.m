function [p0,d,res]=svdfit(D,extraP)
% fit dots in a line in 3D
% D=[x_clustered(clustered==i) y_clustered(clustered==i) z_clustered(clustered==i)];
avg = mean(D, 1); %mean of D along the column
subtracted = bsxfun(@minus, D, avg);% elementwise subtraction of avg from D
[U, S, V] = svd(subtracted);
direction = V(:, 1); %the first principal component of the centered data in matrix D using PCA and stores it in the direction vector
p0 = avg;
d = direction;


% residual
res=sum((subtracted*V(:,end)).^2); %fitting mean centered data using last principle component (not used in the code



% if exist('extraP')
%     avg2 = mean([D;extraP]), 1);
%     subtracted2 = bsxfun(@minus, [D;extraP], avg);
%     [U2, S2, V2] = svd(subtracted);
%     sum((subtracted2*V2(:,end)).^2);
% end