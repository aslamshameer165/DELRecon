% Display or remove the center of mass of each cluster from the figure
[az el]=view;
XYZlim=[get(gca,'Xlim');get(gca,'Ylim');get(gca,'Zlim')];
set(gca,'Xlim',XYZlim(1,:));
set(gca,'Ylim',XYZlim(2,:));
set(gca,'Zlim',XYZlim(3,:));
view(az,el);
if strcmp(get(S.GS,'Visible'),'off')
    set(S.GS,'Visible','on');
else
    set(S.GS,'Visible','off');
end
