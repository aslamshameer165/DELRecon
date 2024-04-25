% Add new contact

[xd,~]=ginput(1);% Get the user defined position of the new contact

% Determine the contact position based on the plot
[~,xd]=min(abs(linesp-repmat(xd,length(linesp),1)));
xd=linesp(xd);

% Update the plots and Contact variable
hold on;
GS2=[GS2;repmat(xd,1,3)];
delete(S.pt)

hold on
if Contacts.contNum(pEle)==length(GS2)
Contacts.contNum(pEle)=Contacts.contNum(pEle)+1;
    
    S.pt=scatter(GS2(:,1),zeros(Contacts.contNum(pEle),1)+(Contacts.yrange(pEle,1)+Contacts.yrange(pEle,2))/2,'MarkerEdgeColor',[0 1 0],...
              'MarkerFaceColor',[0 .8 0]); 

else
    Contacts.contNum(pEle)=Contacts.contNum(pEle)+1;
S.pt=scatter(GS2(:,1),zeros(Contacts.contNum(pEle),1)+(Contacts.yrange(pEle,1)+Contacts.yrange(pEle,2))/2,'MarkerEdgeColor',[1 0 0],...
              'MarkerFaceColor',[.8 0 0]);
end
hold off
text_box1.String=num2str(length(GS2));