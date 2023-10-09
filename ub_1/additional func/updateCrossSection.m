
% Updating the latest Real entry using user defined cross section
% properties: A, Izz, Iyy, and Ixx
% Inp.Real(size(Inp.Real,1)-numOfSecElems+1:size(Inp.Real,1),2) = A;
% Inp.Real(size(Inp.Real,1)-numOfSecElems+1:size(Inp.Real,1),3) = Izz;
% Inp.Real(size(Inp.Real,1)-numOfSecElems+1:size(Inp.Real,1),4) = Iyy;
% Inp.Real(size(Inp.Real,1)-numOfSecElems+1:size(Inp.Real,1),9) = Ixx;

% ansEntry has the previous layer  so with update section if senod layer is
% also added
    
Inp.Real(layerIndex(layerID,:),2) = A;
Inp.Real(layerIndex(layerID,:),3) = Izz;
Inp.Real(layerIndex(layerID,:),4) = Iyy;
Inp.Real(layerIndex(layerID,:),9) = Ixx;