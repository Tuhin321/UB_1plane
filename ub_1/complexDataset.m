complexUB(:,1)=metadata(:,1).*exp(1i*metadata(:,2));

% Check abs(metadataComplex)
%Request.UBResp.OutputNodes=[2 9 22 28]; % List of Output Nodes



complexUB(:,2)= metadata(:,4)+1i*metadata(:,5);
complexUB(:,3)= metadata(:,6)+1i*metadata(:,7);
complexUB(:,4)= metadata(:,12)+1i*metadata(:,13);
complexUB(:,5)= metadata(:,14)+1i*metadata(:,15);
realonly=real(complexUB);

% writematrix(complexUB,'Datasetreal_v2_2.csv') 
writematrix(realonly,'Datasetreal_v2_3.csv') 