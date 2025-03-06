clear;

Delta = 2.5e-4;
D = 1;
ozin = [-1;1]*D;
RhoV2in = [1;1]*(Delta/pi);

% NRG parameters
Lambda = 2.5;
N = 55;
Nkeep = 500;

[ff,gg] = CLD_Prac(ozin,RhoV2in,Lambda,N);

% impurity Hamiltonian parameters
U = 4e-3;
epsd = -U/2;
J = 8*Delta*D/(pi*U);

% operators for first bath site (0th site in Wilson chain)
[F,Z,S,I] = getLocalSpace('FermionS','Acharge,SU2spin','NC',1);
F.info.itags = {'L0','L0*','op*'};
Z.info.itags = {'L0','L0*'};
S.info.itags = {'L0','L0*','op*'};
I.E.info.itags = {'L0','L0*'};

% operators for impurity site ('-1'th site in Wilson chain)
S_imp = S;
I_imp = getsub(I.E,2);
S_imp.info.itags = {'imp','imp*','op*'};
I_imp.info.itags = {'imp','imp*'};

% Particle number operator in first bath site (0th site in Wilson chain)
NF = contract(F,'!2*',F,'!2');

% ket tensor for the impurity + first bath site
A0 = getIdentity(I_imp,2,I.E,2,'H0*',[1,3,2]);

% Local Hamiltonian
H0 = 2*J*contract(A0,'!2*',{S,'!2*',{S_imp,A0}});
H0 = H0 + gg(1)*contract(NF,{A0,'!2*',A0,'!3'});

% Iterative diagonalization
Inrg = NRG_IterDiag_Prac(H0,A0,Lambda,ff(2:end),F,gg(2:end),NF,Z,Nkeep);

% Energy flow diagram
Eshow = 3; % energy window to show (from 0 to Eshow)

Eeven = cell(1,floor((numel(Inrg.EK)+1)/2));
Q_num_even = cell(1,floor((numel(Inrg.EK)+1)/2));

for itN = (1:2:numel(Inrg.EK))
    for it = (1:numel(Inrg.EK{itN}.data))
        Eeven{(itN+1)/2} = cat(1 , Eeven{(itN+1)/2} , Inrg.EK{itN}.data{it}.');
        Q_num_even{floor((itN+1)/2)} = cat(1 , Q_num_even{floor((itN+1)/2)} , cat(2 , Inrg.EK{itN}.Q{1}(it,1)*ones(numel(Inrg.EK{itN}.data{it}),1) , ...
            Inrg.EK{itN}.Q{1}(it,2)*ones(numel(Inrg.EK{itN}.data{it}),1) , Inrg.EK{itN}.data{it}.'));
    end
    [Eeven{(itN+1)/2}, ids] = sort(Eeven{(itN+1)/2},'ascend');
    Q_num_even{floor((itN+1)/2)} = Q_num_even{floor((itN+1)/2)}(ids,:);

    Q_num_even{floor((itN+1)/2)} = Q_num_even{floor((itN+1)/2)}( Eeven{floor((itN+1)/2)}(:,1) <= Eshow, :);
end


%since we start from A0; consider the step for H0 as 0, i.e., even
%Eeven = Inrg.EK_arr(1:2:end);
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
Eodd = cell(1,floor(numel(Inrg.EK)/2));
Q_num_odd = cell(1,floor(numel(Inrg.EK)/2));

for itN = (2:2:numel(Inrg.EK))
    for it = (1:numel(Inrg.EK{itN}.data))
        Eodd{floor(itN/2)} = cat(1 , Eodd{itN/2} , Inrg.EK{itN}.data{it}.');
        Q_num_odd{floor(itN/2)} = cat(1 , Q_num_odd{floor(itN/2)} , cat(2 , Inrg.EK{itN}.Q{1}(it,1)*ones(numel(Inrg.EK{itN}.data{it}),1) , ...
            Inrg.EK{itN}.Q{1}(it,2)*ones(numel(Inrg.EK{itN}.data{it}),1) , Inrg.EK{itN}.data{it}.'));
    end
    [Eodd{itN/2}, ids] = sort(Eodd{itN/2},'ascend');
    Q_num_odd{floor(itN/2)} = Q_num_odd{floor(itN/2)}(ids,:);

    Q_num_odd{floor(itN/2)} = Q_num_odd{itN/2}( Eodd{floor(itN/2)}(:,1) <= Eshow , :);
end

%Eodd = Inrg.EK_arr(2:2:end);
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