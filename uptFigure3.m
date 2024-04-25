% Display or remove the tip of electrode handles from the figure
[az el]=view;
XYZlim=[get(gca,'Xlim');get(gca,'Ylim');get(gca,'Zlim')];
set(gca,'Xlim',XYZlim(1,:));
set(gca,'Ylim',XYZlim(2,:));
set(gca,'Zlim',XYZlim(3,:));
view(az,el);
if strcmp(get(S.seed,'Visible'),'off')
    set(S.seed,'Visible','on');
else
    set(S.seed,'Visible','off');
end
