function h=plotBrainandElec(ima_mask_erode2,ima_clustered,x_clustered,y_clustered,z_clustered,x,y,z,clustered,GS_clustered,clr,showCod)

h.fig = figure('DeleteFcn','doc datacursormode','Position',[100 100 1000 500]);
h.ax1 = subplot(1, 2, 1); hold on; % Subplot to display the 3D brain with the initial clusters
bb=isosurface(permute(ima_mask_erode2,[2,1,3])); %to extract a 3D isosurface from cortical surface strip
h.p=patch(bb); % Create patch graphics object to display isosurface
% Define the face and edge colors for the isosurface
h.p.FaceColor=[.5 .5 .5];
h.p.EdgeColor=[.7 .7 .7];
reducepatch(h.p,4000);% Reduce the number of faces of patch object to 4000
% Plot the initial clusters
h.ax1handleplot=scatter3(x_clustered, y_clustered,z_clustered,36,clr(clustered,:),'Marker','.');
h.ax1datplot=scatter3(x,y, z,36,'Marker','.','MarkerEdgeColor',[0 0 1]);
view(3);axis vis3d, box on;rotate3d on;
xlabel('x'),ylabel('y'),zlabel('z');

% Plot only initial clusters in another subplot
h.ax2 = subplot(1,2,2);
h.ax2datplot=scatter3(x,y, z,36,'Marker','.','MarkerEdgeColor',[0 0 1]);
hold on;
h.ax2handleplot=scatter3(x_clustered, y_clustered, z_clustered,36,clr(clustered,:),'Marker','.');
h.ax2GMplot=scatter3(GS_clustered(:,1),GS_clustered(:,2),GS_clustered(:,3),100,clr(unique(clustered),:),'Marker','o','LineWidth',3);
% view(3);axis vis3d, box off;rotate3d on;
% xlabel('x'),ylabel('y'),zlabel('z');

h.ax2.XLim=h.ax1.XLim;
h.ax2.YLim=h.ax1.YLim;

% Link properties of two axes
Link = linkprop([h.ax1, h.ax2], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);

% Set up the properties of data cursor object using custom dataTipFormat function
dcm_obj = datacursormode(h.fig);
set(dcm_obj,'UpdateFcn',{@dataTipFormat,ima_clustered,showCod});

% Create interactive data cursor annotation for the data points
for i=1:max(unique(clustered))
    hdtip(i)= dcm_obj.createDatatip(h.ax2handleplot);
    %     set(hdtip(i), 'MarkerSize',5, 'MarkerFaceColor','none', 'MarkerEdgeColor','r', 'Marker','o', 'HitTest','off');
    curCL = find(clustered==i);
    set(hdtip(i), 'Position', [x_clustered(curCL(1)) y_clustered(curCL(1)) z_clustered(curCL(1))]);
    updateDataCursors(dcm_obj);
    pause(0.2);
end