% Lower electrode length range
[xr1,~]=ginput(1);% get the user defined lower range

% Determine the range value based on the plot
[~,xr1]=min(abs(linesp-repmat(xr1,length(linesp),1)));
xr1=linesp(xr1);

% Update the Contacts variable with new range value
Contacts.xrange(pEle,1)=xr1;

% Reclustering to get individual contact
clear indx xx GS2 clusters2 linespGS
    indx=linedist<2 & linesp<=Contacts.xrange(pEle,2) & linesp>=Contacts.xrange(pEle,1);
    
    xx=sub2ind(size(ima),P(indx,1),P(indx,2),P(indx,3));
    [GS2,clusters2]=  clustering(Contacts.contNum(pEle),linesp(indx),ima(xx));
    
    % Update the plots
    linespGS = svd3D_1D(p0(pEle,:), d(:,pEle)',GS_seed(pEle,:)) ;
    hold on
    delete(S.sc)
    S.sc=scatter(linesp(linedist<20),bbb(linedist<20)-uu,1,repmat(linedist(linedist<20)./max(linedist),1,3));
    ylim([uu,max(bbb(linedist<20))-uu])
    delete(S.pt)
    S.pt=scatter(GS2(:,1),zeros(Contacts.contNum(pEle),1)+(Contacts.yrange(pEle,1)+Contacts.yrange(pEle,2))/2,'r');
delete(S.ln1)
delete(S.ln2)
delete(S.ln1y)
delete(S.ln2y)
    S.ln1=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,1)],[0 max(bbb(linedist<20))],'r');
    S.ln2=plot([Contacts.xrange(pEle,2) Contacts.xrange(pEle,2)],[0 max(bbb(linedist<20))],'r');
    S.ln1y=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,2)],[Contacts.yrange(pEle,1) Contacts.yrange(pEle,1)],'b');
    S.ln2y=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,2)],[Contacts.yrange(pEle,2) Contacts.yrange(pEle,2)],'b');
hold off;