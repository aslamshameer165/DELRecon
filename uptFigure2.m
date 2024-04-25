% To display or remove voxels beyond the cortical surface strip from the figure
[az el]=view;
XYZlim=[get(gca,'Xlim');get(gca,'Ylim');get(gca,'Zlim')];
set(gca,'Xlim',XYZlim(1,:));
set(gca,'Ylim',XYZlim(2,:));
set(gca,'Zlim',XYZlim(3,:));
view(az,el);
if strcmp(get(S.all,'Visible'),'off')
    set(S.all,'Visible','on');
else
    set(S.all,'Visible','off');
end
