% Zoom out the plots by 100 units
% update all the plots   
uu=uu-100;
    hold on
   delete(S.sc)
   delete(S.pt)
   delete(S.ln1)
   delete(S.ln2)
   delete(S.ln1y)
   delete(S.ln2y)
    S.sc=scatter(linesp(linedist<20),bbb(linedist<20)-uu,1,repmat(linedist(linedist<20)./max(linedist),1,3));
    ylim([uu,max(bbb(linedist<20))-uu])
    S.pt=scatter(GS2(:,1),zeros(Contacts.contNum(pEle),1)+(Contacts.yrange(pEle,1)+Contacts.yrange(pEle,2))/2,'r');
    S.ln1=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,1)],[0 max(bbb(linedist<20))],'r');
    S.ln2=plot([Contacts.xrange(pEle,2) Contacts.xrange(pEle,2)],[0 max(bbb(linedist<20))],'r');
    S.ln1y=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,2)],[Contacts.yrange(pEle,1) Contacts.yrange(pEle,1)],'b');
    S.ln2y=plot([Contacts.xrange(pEle,1) Contacts.xrange(pEle,2)],[Contacts.yrange(pEle,2) Contacts.yrange(pEle,2)],'b');
    hold off