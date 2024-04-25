% Reset the figure 
delete(S.sc)
   uu=0;
    
   % Plot all the data points
    S.sc=scatter(linesp(linedist<20),bbb(linedist<20)-uu,1,repmat(linedist(linedist<20)./max(linedist),1,3));
    ylim([uu,max(bbb(linedist<20))-uu])
    % Reset all the range values to the default
        Contacts.xrange(pEle,1)=min(linesp);
        Contacts.xrange(pEle,2)=max(linesp);
        Contacts.yrange(pEle,1)=min(bbb(linedist<20));
        Contacts.yrange(pEle,2)=max(bbb(linedist<20));
        
   % Recluster to get individual contact labels
   indx=linedist<2 & linesp<=Contacts.xrange(pEle,2) & linesp>=Contacts.xrange(pEle,1) & bbb(linedist<20)<=Contacts.yrange(pEle,2) & bbb(linedist<20)>=Contacts.yrange(pEle,1);

   xx=sub2ind(size(ima),P(indx,1),P(indx,2),P(indx,3));
    [GS2,clusters2]=  clustering(Contacts.contNum(pEle),linesp(indx),ima(xx));

    % Update all the plots
    linespGS = svd3D_1D(p0(pEle,:), d(:,pEle)',GS_seed(pEle,:)) ;
    delete(S.pt)
    delete(S.ln1)
    delete(S.ln2)
    delete(S.ln1y)
    delete(S.ln2y)
    hold on
    S.pt=scatter(GS2(:,1),zeros(Contacts.contNum(pEle),1)+(Contacts.yrange(pEle,1)+Contacts.yrange(pEle,2))/2,'r');
    S.ln1=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,1)],[0 max(bbb(linedist<20))],'r');
    S.ln2=plot([Contacts.xrange(pEle,2) Contacts.xrange(pEle,2)],[0 max(bbb(linedist<20))],'r');
    S.ln1y=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,2)],[Contacts.yrange(pEle,1) Contacts.yrange(pEle,1)],'b');
    S.ln2y=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,2)],[Contacts.yrange(pEle,2) Contacts.yrange(pEle,2)],'b');
    hold off