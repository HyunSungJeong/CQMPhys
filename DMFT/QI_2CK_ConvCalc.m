function [ocont, DynSusc, ThermExp] = QI_2CK_ConvCalc(DMFTdata,Nkeep,Lambda,nz,U,t_0,phi_div_pi,T,varargin)
% <Description>
% 
% computes physical quantities (spectral functions, dynamic sysceptibilities, etc.)
% from the converged hybridization function from DMFT calculations (obtained by running 'DMFT_QI_2CK_lat.m')
%
% <Input>
% DMFTdata : [char] path to the folder where DMFT results are stored
%           (e.g. '/data/hyunsung/DMFT_QI_2CK_lat/U=0_V=0.5_t_0=0.5_phi_div_pi=0.25_T=1e-10_ndfix=5_Nkeep=3000')
% Nkeep : [integer] number of kept states in each NRG itereation
% Lambda : [numeric] logarithmic discretization parameter
% nz : [integer] number of z-shifts
% U : [numeric] on-site Coulomb repulsion
% t_0 : [numeric] intra-cell hopping amplitude t_0
% phi_div_pi : [numeric] distortion angle phi/pi
% T : [numeric] temperature
%
% <Option>
% 'DynSusc', Ops1, Ops2, Z, zflag, cflag, g : 
%           Calculate dynamic susceptibilities
%           Ops1, Ops2 : [(array of) QSpace] Ops1 and Ops2 are local operators to be fed to 'getAdisc' for dynamic susceptibility calculations.
%                               The requirements for Ops1 and Ops2 are identical to that specified in the documentation of 'getAdisc'
%           Z : [(array of) QSpace] Fermionic sign operator at site s00.
%           zflag : [numeric (vector)] To use (1) or not to use (0) fermionic sign change for each correlator.
%           cflag : [numeric (vector)] Sign factor(s) for each commutator.
%                           [Op1(n), Op2(n)']_{+ or -}. +1 for fermionic (= anti-commuting) operators, 
%                           -1 for bosonic (= commuting) operators.
%           g : [numeric] Parameter for secondary linear broadening kernel. (\gamma in Lee2016.)
%
% 'useDMFTchain', .. : [Logical] if true, use Wilson chain parameters used in the final iteration of DMFT calculation
%                               if false, use Wilson chain parameters obtained by discretizing the converged hybridization function
% 'doZLD', ... : [cell array] Options for 'doZLD' (e.g. {'Nfit',30})
%                           (Default: {'Nfit', round(-2*log(1e-12)/log(Lambda))})
% 'NRG_SL', .. : [cell array] Options for 'NRG_SL' (e.g. {'deps',1e-10})
%                           'Lambda' is already set. 
%                           (Default: {'Lambda',Lambda,'deps',1e-10})                        
% 'getAdisc', .. : [cell array] Options for 'getAdisc' (e.g. {'emin',1e-10})
% 'getAcont', .. : [cell array] Options for 'getAcont' (e.g. {'alphaz',1,'emin',1e-10})
% 'ThermExp', .. : [(array of) QSpace] Local operators for thermal expectation value calculations using fdmNRG
%               They must act either on the local physical spaces '[Ls]00[*]', on the kept space 'K00[*]', or on the discarded space 'D00[*]'
%
% <Output>
% 'DynSusc' : [cell array] Dynamic susceptibilities corresponding to Ops1 and Ops2
% 'ThermExp' : [numeric array] Thermal expectation values of the input operators in the fdmNRG framework

%% parse inputs and options

% input sanity check
if ~iesqual(class(DMFTdata),'char')
    error('DMFTdata must be a path to the DMFT results');
end

if ~isnumeric(U)
    error('U must be a real number');
end

if ~isnumeric(t_0)
    error('t_0 must be a real number');
end

if ~isnumeric(phi_div_pi)
    error('phi/pi must be a real number');
end

while ~isempty(varargin)
    switch varargin{1}
        case 'DynSusc'
            if isequal(class(varargin{2}),'QSpace')
                Ops1 = varargin{2};
            else
                error('Ops1 must be a QSpace object');
            end

            if isequal(class(varargin{3}),'QSpace')
                Ops2 = varargin{3};
            else
                error('Ops2 must be a QSpace object');
            end
            
            if isequal(class(varargin{4}),'QSpace')
                Z = varargin{4};
            else
                error('Z must be a QSpace object');
            end

            if ~isnumeric(varargin{5})
                error('zflag must be a numeric vector');
            elseif ~isnumeric(varargin{6})
                error('cflag must be a numeric vector');
            elseif numel(varargin{5}) ~= numel(varargin{6})
                error('zflag and cflag must have same lengths');
            else
                zflag = varargin{5};
                cflag = varargin{6};
            end

            if isnumeric(varargin{7})
                g = varargin{7};
            else
                error('g must be a positive real number');
            end

            varargin(1:7) = [];

        case 'useDMFTchain'
            if isequal(class(varargin{2}), 'logical')
                useDMFTchain = varargin{2};
            else
                error('useDMFTchain must be either true or false');
            end
            varargin(1:2) = [];

        case 'doZLD'
            if isequal(class(varargin{2}),'cell')
                opt_doZLD = varargin{2};
            else
                error('doZLD must be a cell array containing the options for the function doZLD');
            end

        case 'NRG_SL'
            if isequal(class(varargin{2}),'cell')
                opt_NRG_SL = varargin{2};
            else
                error('NRG_SL must be a cell array containing the options for the function NRG_SL');
            end

        case 'getAdisc'
            if isequal(class(varargin{2}),'cell')
                opt_getAdisc = varargin{2};
            else
                error('getAdisc must be a cell array containing the options for the function getAdisc');
            end

        case 'getAcont'
            if isequal(class(varargin{2}),'cell')
                opt_getAcont = varargin{2};
            else
                error('getAcont must be a cell array containing the options for the function getAcont');
            end

        case 'ThermExp'
            if isequal(class(varargin{2}),'QSpace')
                Ops_ThermExp = varargin{2};
            else
                error('ThermExp (array of) QSpace object(s)');
            end

        otherwise
            if ischar(varargin{1})
                error(['ERR: Unknown input ',varargin{1}]);
            else
                disp2(varargin{1});
                error('ERR: Unknown input');
            end
    end
end

    %% system parameters
    phi = pi*phi_div_pi;    % Distortion angle \phi
    t = t_0*cos(phi);       % intra-cell hopping amplitude t
    t_p = t_0*sin(phi);     % intra-cell hopping amplitude t'

    %% define local operators

    % local operators in the even/odd basis
    [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge,SU2spin', 'NC', 3);
    [FF_L_eo, ZF_L_eo, EF_L] = setItag('L00', 'op', FF, ZF, IF.E);   % local operators in the L00(impurity) site. [2e, 3, 2o] orbitals
        
    [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge,SU2spin', 'NC', 2);
    [FF_s_eo, ZF_s_eo, EF_s] = setItag('s00', 'op', FF, ZF, IF.E);   % local operators in the s00(first bath) site. [1e, 1o] oribitals

    % number operators in the even/odd basis
    NF_L_eo = quadOp(FF_L_eo, FF_L_eo, []);    % number operators in the L00(impurity) site
    NF_s_eo = quadOp(FF_s_eo, FF_s_eo, []);    % number opators in the s00(first bath) site

    % annihilation ops. in the 'original' basis
    FF_L_orig = QSpace(1,3);        % annihilation ops. for orbital # 2,3,4 (defined at the L00 site)
    FF_s_orig = QSpace(1,2);        % annihilation ops. for oribtal # 1,5 (defined at the s00 site)
        
    FF_L_orig(1) = (FF_L_eo(1) + FF_L_eo(3)) / sqrt(2);
    FF_L_orig(2) = FF_L_eo(2);
    FF_L_orig(3) = (FF_L_eo(1) - FF_L_eo(3)) / sqrt(2);
    FF_s_orig(1) = (FF_s_eo(1) + FF_s_eo(2)) / sqrt(2);
    FF_s_orig(2) = (FF_s_eo(1) - FF_s_eo(2)) / sqrt(2);

    % number operators in the 'original' basis
    NF_L_orig = quadOp(FF_L_orig, FF_L_orig, []);       % number operators in the L00(impurity) site
    NF_s_orig = quadOp(FF_s_orig, FF_s_orig, []);       % number opators in the s00(first bath) site

    %% define impurity Hamiltonian

    % isometry connecting impurity and first bath site
    A0 = getIdentity(EF_L, 2, EF_s, 2, 'K00*', [1 3 2]);

    % Hamiltonian - interacting part(HU)
    HU_L = QSpace();    % part of HU that is defined on the L00 site
    for itL = 1:3
        HU_L = HU_L + (U/2)*(contract(NF_L_orig(itL),'!1',NF_L_orig(itL)) - NF_L_orig(itL));
    end
    HU_L = HU_L - (U/2)*sum(NF_L_orig);

    HU_s = QSpace();    % part of HU that is defined on the s00 site
    for its = 1:2
        HU_s = HU_s + (U/2)*(contract(NF_s_orig(its),'!1',NF_s_orig(its)) - NF_s_orig(its));
    end
    HU_s = HU_s - (U/2)*sum(NF_s_orig);

    HU = contract(A0,'!2*',{HU_s,A0}) + contract(A0,'!2*',{HU_L,A0});   % full interacting part Hamiltonian
    HU = HU + 1e-40*contract(A0,'!2*',A0);                              % add infinitesimal*identity to avoid empty qspace when U=0

    % Hamiltonian - intra-cell hopping(H_hop)
    H_hop = sqrt(2)*t_0*contract(FF_L_eo(1),'!2*',FF_L_eo(2),'!2');
    H_hop = contract(A0,'!2*',{H_hop,A0});
    H_hop = H_hop + (t + t_p)*contract(A0,'!2*',{FF_s_eo(1),{FF_L_eo(1),'!2*',A0}});
    H_hop = H_hop + (t - t_p)*contract(A0,'!2*',{FF_s_eo(2),{FF_L_eo(3),'!2*',A0}});
    H_hop = H_hop + H_hop';

    H0 = HU + Hhop;     % full impurity Hamiltonain

    %% NRG calculation from converged hybridization function
    
    if useDMFTchain     % use Wilson chain parameters from the last DMFT iteration
        tmp = load([DMFTdata,'/ffs.mat']);
        field = getfield(tmp);
        ffs = getfield(tmp,field{1});
        ff = ffs{end};

        tmp = load([DMFTdata,'/ggs.mat']);
        field = getfield(tmp);
        ggs = getfield(tmp,field{1});
        gg = ggs{end};

        tmp = load([DMFTdata,'/dffs.mat']);
        field = getfield(tmp);
        dffs = getfield(tmp,field{1});
        dff = dffs{end};

        tmp = load([DMFTdata,'/dggs.mat']);
        field = getfield(tmp);
        dggs = getfield(tmp,field{1});
        dgg = dggs{end};

        nz = size(ff,2);

    else        % obtain Wilson chain parameters by discretizing the converged hybridization function
        tmp = load([DMFTdata,'/RhoV2in.mat']);
        field = fieldnames(tmp);
        RV2all = getfield(tmp,field{1});
        RV2conv = RV2all(:,:,end);

        tmp = load([DMFTdata,'/ocont.mat']);
        field = fieldnames(tmp);
        ocont = getfield(tmp,field{1});

        Nfit = round(-2*log(1e-12)/log(Lambda));
        [ff, gg, dff, dgg] = doZLD(ocont, RV2conv, Lambda, N, nz, opt_doZLD{:});
    end

    nrgdata = cell(1,nz);

    for itz = (1:nz)
        nrgdata{itz} = NRG_SL([],H0,A0,Lambda,ff{itz}(2:end),Fs,Zs,opt_NRG_SL{:});
        nrgdata{itz} = getRhoFDM(nrgdata{itz},T,'-v','Rdiag',true);   % calculating the full density matrix(FDM)
    end

    %% Calculate physical quantities

    if Calc_DynSusc      % calculate dynamic susceptibilities

        Adiscs = cell(numel(Ops1),nz);      % discrete data
        DynSusc = cell(1,numel(Ops1));    % continuous (i.e., broadened) spectral function
        
        for itz = (1:nz)
            [odisc,Adiscs,sigmak] = getAdisc(nrgdata{itz},Ops1,Ops2,Z,'zflag',zflag,'cflag',cflag,opt_getAdisc{:});
        end
        
        for ita = (1:size(Adiscs,1))
            Adisc = mean(cell2mat(reshape(Adiscs(ita,:),[1 1 nz])),3);
        
            [ocont, DynSusc{ita}] = getAcont(odisc,Adisc,sigmak,g,opt_getAcont{:});
        end
    end

    if Calc_ThermExp

        ThermExp = zeros(1,numel(Ops_ThermExp));

        for it = 1:numel(Ops_ThermExp)
            ThermExp(it) = getEpVal(nrgdata, Ops_ThermExp(it));
        end
    end

end