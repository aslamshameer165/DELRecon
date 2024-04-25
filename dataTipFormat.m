%% Format the data tip updates
function txt = dataTipFormat(~,event_obj,ima_clustered,showCod)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
if showCod==1
    txt = {['Cor: ',[num2str(floor(pos(1))) ', ' num2str(floor(pos(2))) ', ' num2str(floor(pos(3)))]],...
        ['Cl-',num2str(ima_clustered(round(pos(1)),round(pos(2)),round(pos(3))))]};
else
    txt = {['Cl-',num2str(ima_clustered(round(pos(1)),round(pos(2)),round(pos(3))))]};
end