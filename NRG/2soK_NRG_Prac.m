clear;

Delta = 2.5e-4;
D = 1;
ozin = [-1;1]*D;
RhoV2in = [1;1]*(Delta/pi);

% NRG parameters
Lambda = 2.5;
N = 55;
Nkeep = 300;

[ff,gg] = CLD_Prac(ozin,RhoV2in,Lambda,N);

% Impurity Hamiltonian parameters
J0 = 0.1;
I0 = 0.01;
K0 = 0.1;

% Operators for 0th wilson chain site
[F,Z,J_sp,I] = getLocalSpace('FermionS','Acharge,SU2spin,SU2channel','NC',2);
J_sp = J_sp/2;
F.info.itags = {'L00','L00*','op*'};
Z.info.itags = {'L00','L00*'};
J_sp.info.itags = {'L00','L00*','op*'};
I.E.info.itags = {'L00','L00*'};

F1 = getsub(F,1);   
S1 = contract(F1,'!2*',F1);     % [ -1 1 1 ; -1 1 1 ]
S1.data{1} = J_sp.data{1};

F2 = getsub(F,2);
S2 = contract(F2,'!2*',F2);     % [ 0 0 2 ; 0 0 2 ]
S2.data{1} = J_sp.data{2};

F3 = getsub(F,[4,5]);
S3 = contract(F3,'!2*',F3);
S3.data{1} = J_sp.data{3};      % [ 1 1 1 ; 1 1 1 ]

J_orb = S1 + S2 + S3;
J_orb.data = S_sp.data;

F_imp = getsub(F,1);
Z_imp = getsub(Z,2);
I_imp = getsub(I.E,2);
F_imp.info.itags = {'imp','imp*','op*'};
Z_imp.info.itags = {'imp','imp*'};
I_imp.info.itags = {'imp','imp*'};

A = getIdentity(getsub(F,[2,3,5]),2,getsub(F,[2,3,5]),3);

% Particle number operator in first bath site (0th site in Wilson chain)
NF = contract(F,'!2*',F,'!2');

% ket tensor for the impurity + first bath site
A0 = getIdentity(I_imp,2,I.E,2,'H0*',[1,3,2]);

% Local Hamiltonian
H0 = J0*contract(A0,'!2*',{J_sp,'!2*',{S_sp,A0}});
H0 = H0 + K0*contract(A0,'!2*',{J_orb,'!2*',{S_orb,A0}});
H0 = H0 + I0*contract()
H0 = H0 + gg(1)*contract(NF,{A0,'!2*',A0,'!3'});

% Iterative diagonalization
Inrg = NRG_IterDiag_Prac(H0,A0,Lambda,ff(2:end),F,gg(2:end),NF,Z,Nkeep);

% Energy flow diagram
Eshow = 3; % energy window to show (from 0 to Eshow)

%since we start from A0; consider the step for H0 as 0, i.e., even
Eeven = Inrg.EK_arr(1:2:end);
Eeven = cellfun(@(x) x(x <= Eshow), Eeven, 'UniformOutput',0); 
% 1D cell array, 1D cell array, each cell is vector containing kept energies not exceeding Eshow for each even step
maxEeven = max(cellfun('prodofsize',Eeven));    % maximum dimension of vectors in Eeven
Eeven = cellfun(@(x) [x;nan(maxEeven-numel(x),1)], Eeven, 'UniformOutput',0);
% 1D cell array same as above, but each cells are now equal size
% nan filled so each cells are vectors of equal size maxEeven
% Can now be concatenated into a single matrix
Eeven = cell2mat(Eeven).';  
% Eeven concatenated into a single matrix
% each row vector containing energies below Eshow 
% for corresponding even iteration step

% same for odd iterations
Eodd = Inrg.EK_arr(2:2:end);
Eodd = cellfun(@(x) x(x <= Eshow), Eodd, 'UniformOutput',0);
maxEodd = max(cellfun('prodofsize',Eodd));
Eodd = cellfun(@(x) [x;nan(maxEodd-numel(x),1)], Eodd, 'UniformOutput',0);
Eodd = cell2mat(Eodd).';

figure;
% upper panel
subplot(2,1,1);
plot((1:2:numel(Inrg.EK))-1,Eeven,'LineWidth',1);
xlabel('Even Iterations');
xlim([0 numel(Inrg.EK)-1]);
ylim([0 Eshow]);
set(gca,'LineWidth',1,'FontSize',13);

%lower panel
subplot(2,1,2);
plot((2:2:numel(Inrg.EK))-1,Eodd,'LineWidth',1);
xlabel('Odd Iterations');
xlim([0 numel(Inrg.EK)-1]);
ylim([0 Eshow]);
set(gca,'LineWidth',1,'FontSize',13);