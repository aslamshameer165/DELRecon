function xyz = svd1D_3d(p0,d,linesp)

xyz = linesp*d + repmat(p0,length(linesp),1);
