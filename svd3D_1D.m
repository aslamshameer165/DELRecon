function linesp = svd3D_1D(p0,d,xyz)

linesp = (xyz-repmat(p0,size(xyz,1),1))/d; % The unit vector from SVD center point p0 to the xyz along d

