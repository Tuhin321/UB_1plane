
% Storing the ans entry
ansEntry = ans;

% colors for 8 materials
newcolors=[0 1 1; 0.5 0.5 0.5; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980; 0.3010 0.7450 0.9330; 0 0.4470 0.7410;0.4940 0.1840 0.5560; 0.6350 0.0780 0.1840];
% newcolors=[0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410;0.4940 0.1840 0.5560; 0.6350 0.0780 0.1840];
% Section length
Lsection = ansEntry(1);


% Max diameter of the section
Dmax = ansEntry(end-1);


% Automatic section disctretization
numOfSecElems = round(2*Lsection/Dmax);
if numOfSecElems < 1 % minimum of one element per section length
    numOfSecElems = 1;
end


% Number of element layers
numOfElemLayers = (length(ansEntry)-2)/2;


% Adding nodes
newNodes = ((Inp.Node(end,1)+1):(Inp.Node(end,1)+numOfSecElems))';
Xcoords = linspace(Inp.Node(end,2),Inp.Node(end,2)+ansEntry(1),numOfSecElems+1)';

Inp.Node = [Inp.Node; [newNodes Xcoords(2:end) ones(numOfSecElems,1)*Inp.Node(end,3) ones(numOfSecElems,1)*Inp.Node(end,4)]];


% Adding real constants and elements
if ~isfield(Inp,'Real')
    Inp.Real = [];
end

if ~isfield(Inp,'Elem')
    Inp.Elem = [];
end

layerIndex = [];

for layer = 1:numOfElemLayers
    DoElement = ansEntry(2+2*(layer-1)+1);
    if layer == 1
        DiElement = ansEntry(2+2*(layer-1));
    else
        DiElement = ansEntry(2+2*(layer-1)-1);
    end
    MatID = ansEntry(2+2*(layer-1)+2);
    nuElement = Inp.Mat(MatID,3);
    for element = 1:numOfSecElems
        Inp.Real = [Inp.Real; size(Inp.Real,1)+1 calculateRealConstants(DoElement,DiElement,nuElement)];
        Inp.Elem = [Inp.Elem; size(Inp.Elem,1)+1 newNodes(element)-1 newNodes(element) MatID size(Inp.Real,1)];
        Inp.ElemColor(size(Inp.Elem,1),:)= newcolors(MatID,:);
        layerIndex(layer,element) = size(Inp.Real,1);

    end
end
