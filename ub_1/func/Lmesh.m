function [Node,Elem,MaxNodeNro,MaxElemNro]=Lmesh(startp,endp,Nelem,MatID,RealID,nodestart,elemstart,start_node)
% This function generates Node and Elem matrices for a line
% [Node,Elem,MaxNodeNro,MaxElemNro]=Lmesh(startp,endp,Nelem,MatID,RealID,nodestart,elemstart,start_node)
% 
% Output:
% -------
% Node = Node matrix
% Elem = Element matrix
% MaxNodeNro = maximum node number
% MaxElemNro = maximum element number

% Input:
% ------
% startp = column vector of line start point coordinates
% endp   = column vector of line end point coordinates
%
% Nelem  = number of elements on line
% MatID  = material property number of the elements
% RealID = real constant number of the elements
%
% nodestart = node numbering starting number (default=1)
% elemstart = element numbering starting number (default=1)
% startnode = if first node of the line exists use 0, otherwise 1 (when 0 the node
%             of the first element is not generated)
% 
% Written by Jussi Sopanen LUT/IMVe 2004, jsopanen@lut.fi

%-------------------------------------------------------------------------------------------
% if starting numbering is not given, default is 1
if nargin==6
    nodestart=1;
    elemstart=1;
    start_node=1;
end

% element length
L=norm(endp-startp)/Nelem;
% line length
LL=norm(endp-startp);

% preallocating matrices
if start_node==1; 
    Node=zeros(Nelem+1,4); 
else
    Node=zeros(Nelem,4); 
end
Elem=zeros(Nelem,5);

% preallocating node and element numbers
Innum=nodestart-1;
enum=elemstart-1; 

for ii=1:1:Nelem
    
    % calculating node positions
	InodeLoc=startp+(endp-startp)*(ii-1)*L/LL;
	JnodeLoc=startp+(endp-startp)*(ii)*L/LL;

    % node numbers
    Innum=Innum+1;
    Jnnum=Innum+1;

    % Node matrix
    % If first element in line I node does not exist
    if ii==1 & start_node==1
        Node(ii,:)=[Innum InodeLoc' ];
        Node(ii+1,:)=[Jnnum JnodeLoc' ];
    elseif start_node==0 % If I node exists
        Node(ii+1-1,:)=[Jnnum JnodeLoc' ];
    else % if start_node==1 and ii > 1
        Node(ii+1,:)=[Jnnum JnodeLoc' ];
    end
    
    % element number
    enum=enum+1;
    % Elem matrix
    Elem(ii,:)=[enum Innum Jnnum MatID RealID];
end

% Save maximum node and element number
MaxNodeNro=Jnnum;
MaxElemNro=enum;

