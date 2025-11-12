%% Define parameters

Nkeep = 200;
Lambda = 2.5;
J0 = 0.3;
K_perp = 0.1;
K_z = 0.1;
I_perp = 0;
I_z = 0;
T = 1e-11;

ozin = [-1;1];
RhoV2in = [1;1];
N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);

%% Define local operators

% local operators
[FF,ZF,J_sp,IF] = getLocalSpace('FermionS','Acharge(:),SU2spin','NC',2);
[Fs,Zs,J_sp,Es] = setItag('s00','op',FF(:),ZF,J_sp(:),IF.E);

% J_sp: bath spin op, Tr(S_a S_b) = 1/2 \delta_{ab}
% bath orbital pseudospin op, Tr(T_a T_b) = 1/2 \delta_{ab}

F1F2all = quadOp(Fs(1),Fs(2),'*');
J_orb_plus = -sqrt(2)*getsub(F1F2all, find(all(F1F2all.Q{3}(:,3) == 0, 2)));
F2F1all = quadOp(Fs(2),Fs(1),'*');
J_orb_minus = -sqrt(2)*getsub(F2F1all, find(all(F2F1all.Q{3}(:,3) == 0, 2)));
J_orb_z = (1/2)*( contract(Fs(1),'!2*',Fs(1),'!2') - contract(Fs(2),'!2*',Fs(2),'!2') );

% bath spin-orbital op
J_sporb_plus = (1/sqrt(2))*getsub(F1F2all, find(all(F1F2all.Q{3}(:,3) == 2, 2)));
J_sporb_minus = (1/sqrt(2))*getsub(F2F1all, find(all(F2F1all.Q{3}(:,3) == 2, 2)));
F1F1all = quadOp(Fs(1),Fs(1),'*');
F2F2all = quadOp(Fs(2),Fs(2),'*');
J_sporb_z = 0.5*getsub(F1F1all, find(all(F1F1all.Q{3}(:,3) == 2, 2)));
J_sporb_z = J_sporb_z - 0.5*getsub(F2F2all, find(all(F2F2all.Q{3}(:,3) == 2, 2)));
J_sporb_z = (1/sqrt(2))*J_sporb_z;

Z_imp = getsub(Zs, find(all(Zs.Q{1}(:,1)+Zs.Q{1}(:,2) == -1, 2)));
E_imp = getsub(Es, find(all(Es.Q{1}(:,1)+Es.Q{1}(:,2) == -1, 2)));
S_sp = getsub(J_sp, find(all(J_sp.Q{1}(:,1)+J_sp.Q{1}(:,2) == -1, 2)));

S_orb_plus = getsub(J_orb_plus, find(all(J_orb_plus.Q{1}(:,1)+J_orb_plus.Q{1}(:,2) == -1, 2)));
S_orb_minus = getsub(J_orb_minus, find(all(J_orb_minus.Q{1}(:,1)+J_orb_minus.Q{1}(:,2) == -1, 2)));
S_orb_z = getsub(J_orb_z, find(all(J_orb_z.Q{1}(:,1)+J_orb_z.Q{1}(:,2) == -1, 2)));

S_sporb_plus = getsub(J_sporb_plus, find(all(J_sporb_plus.Q{1}(:,1)+J_sporb_plus.Q{1}(:,2) == -1, 2)));
S_sporb_minus = getsub(J_sporb_minus, find(all(J_sporb_minus.Q{1}(:,1)+J_sporb_minus.Q{1}(:,2) == -1, 2)));
S_sporb_z = getsub(J_sporb_z, find(all(J_sporb_z.Q{1}(:,1)+J_sporb_z.Q{1}(:,2) == -1, 2)));

[Z_imp,E_imp,S_sp,S_orb_plus, S_orb_minus, S_orb_z, S_sporb_plus, S_sporb_minus, S_sporb_z] = ...
    setItag('L00','op',Z_imp(:), E_imp,S_sp(:), S_orb_plus, S_orb_minus, S_orb_z, S_sporb_plus, S_sporb_minus, S_sporb_z);


%% Define the local isometry and Hamiltonian

A0 = getIdentity(E_imp,2,Es,2,'K00*',[1,3,2]);

H0 = J0*contract(A0,'!2*',{J_sp,'!2*',{S_sp,A0}});              % spin-spin
H0 = H0 + (K_perp/2)*contract(A0,'!2*',{J_orb_plus,'!2*',{S_orb_plus,A0}});       % orbital-orbital
H0 = H0 + (K_perp/2)*contract(A0,'!2*',{J_orb_minus,'!2*',{S_orb_minus,A0}});
H0 = H0 + K_z*contract(A0,'!2*',{J_orb_z,'!2*',{S_orb_z,A0}});
A = getIdentity(S_sp,3,S_orb_plus,3,'op*');
H0 = H0 + (I_perp/2)*contract(A0,'!2*',{contract(A,'1,2',contract(S_sp,'!3',S_orb_plus,'!2'),'2,4'),'!2',{A0,J_sporb_plus,'!2*'}});  % spin-orbital
H0 = H0 + (I_perp/2)*contract(A0,'!2*',{contract(conj(A),'1,2',contract(S_sp,'!3*',S_orb_plus,'!1*'),'2,4'),'!2',{A0,J_sporb_plus,'!1'}});
H0 = H0 + I_z*contract(A0, '!2*', contract(S_sp,'!3',S_orb_z,'!2',[1,3,2]), {A0,J_sporb_z,'!2*'});
H0 = H0 + 1e-40*contract(A0,'!2*',A0);


%% Discretization and NRG run

[ff, gg] = doZLD(ozin,RhoV2in,Lambda,N,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));

nrgdata = NRG_SL([],H0,A0,Lambda,ff{1}(2:end),FF,ZF,'Nkeep',Nkeep,'deps',1e-10);
nrgdata = getRhoFDM(nrgdata,T,'-v','Rdiag',true);     % calculating the full density matrix(FDM)

[Etot,Qtot,Qdiff] = plotE(nrgdata,'Emax',10,'legmax',25);       % Data for Eflow diagram
plotE(Etot, Qtot, 'Emax',2,'legmax',15,'Qdiff',[0,0,0]);