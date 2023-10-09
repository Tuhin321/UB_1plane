
% Storing the ans entry
ansEntry = ans;

if ~isfield(Inp,'Bearing')
    Inp.Bearing = [];
end

numOfBrgs = length(Inp.Bearing);

Inp.Bearing(numOfBrgs+1).type = 'Bearing Matrix';
Inp.Bearing(numOfBrgs+1).Inode = size(Inp.Node,1);
Inp.Bearing(numOfBrgs+1).Jnode = ansEntry(3);
Inp.Bearing(numOfBrgs+1).Kb = ansEntry(1)*[0 0 0
                                           0 1 0
                                           0 0 1];
Inp.Bearing(numOfBrgs+1).Cb = ansEntry(2)*[0 0 0
                                           0 1 0
                                           0 0 1];
%
