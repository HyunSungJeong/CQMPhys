clear;
[F,Z,S,I] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1);

%{
% Contraction
F1 = F;
F1.info.itags = {'I','II*','III*'};
NF = contract(F1,'!2*',F1)              % number operator
N2 = contract(F1,'!2*',F1)
N2 = contract(F1,'!1',{F1,'!2*',F1})
N2 = contract(F1,'!2*',{F1,'!1',{F1,'!2*',F1}})
%}

% Generate the isometry that combines the spaces of two legs
%{
E1 = I.E;
E1.info.itags = {'s00','s00*'};
E2 = I.E;
E2.info.itags = {'s01','s01*'};
A = getIdentity(E1,2,E2,2,'A01*',[1,3,2])
%}

% Generate a 1j symbol to invert the direction of the legs
%{
F1 = F;
F1.info.itags = {'s00','s00*','op*'};
F1
I0 = getIdentity(F1,3,'-0')
F1I = contract(F1,'!1',I0,'!2')     % RMT different from the tutorial material!! WHY?
%}

%{
% Automatic trucation of all-zero sectors
M1 = I.E;   M1.info.itags = {'s00','s00*'};
M2 = I.E;   M2.info.itags = {'s01','s01*'};
A = getIdentity(M1,2,M2,2,'A01*',[1,3,2]);
contract(A,'!2*',{M1,'!1',{M2,'!1',A}})     % Identity acting on 16-dim space
A.data{1} = zeros(size(A.data{1}));         % Replacing the first sector RMT data with zero
contract(A,'!2*',{M1,'!1',{M2,'!1',A}})     % Same expression as two lines above, but acts on 15-dim space / all-zero tensor truncated

% If you want the all zero sectors to be kept,
% add the identity operator multiplied by a very small number smaller than double precision
%}


% Eigendecomposition

% Constructing the hopping term acting on two spinful fermionic sites
clear;
[F,Z,S,I] = getLocalSpace('FermionS','Acharge,SU2spin');
% for site s00
F1 = F; F1.info.itags = {'s00','s00*','op*'};
E1 = I.E;   E1.info.itags = {'s00','s00*'};
% for site s01
F2 = F; F2.info.itags = {'s01','s01*','op*'};
E2 = I.E;   E2.info.itags = {'s01','s01*'};
Z2 = Z; Z2.info.itags = {'s01','s01*'};
A = getIdentity(E1,2,E2,2,'A01*',[1,3,2]);

H = contract(A,'!2*',{F1,'!1',{F2,'!2*',{Z2,'!1',A}}}) + ...
    contract(A,'!2*',{F1,'!2*',{Z2,'!1',{F2,'!1',A}}}) + ...
    getIdentity(A,2)*1e-40;     % tight-binding Hamiltonian for two one-orbital spinful sites

celldisp(H.data);
[V,D] = eig(H)
normQS(contract(V,'!2*',V) - getIdentity(A,2))  % checking the unitarity of V
celldisp(D.data)    % D is a QSpace object that has a row vector of eigenvalues as RMT data for each symmetry sector
D2 = diag(D)
celldisp(D2.data)    % D2 = diag(D) is a operator representing a diagonal matrix
[E,Ieig] = eigQS(H,'Nkeep',7)     % using the original MEX function eigQS / Slightly different syntax from the wrap-up eig!
E       % two-column matrix 
        % first column: energy eigenvalues in ascending order
        % second column: multiplet dimensions(degeneracies) for each eigenvlaues
Ieig    % result of MEX functions are MATLAB built-in data type, due to MATLAB policy
% wraping up as QSpace objects(user-defined data type, not a built-in data type)
Ieig.EK = QSpace(Ieig.EK);
Ieig.AK = QSpace(Ieig.AK);
Ieig.EK
Ieig.AK