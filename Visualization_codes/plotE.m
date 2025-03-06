function varargout = plotE (varargin)
% < Description >
% 
% Plot energy flow diagram or entanglement spectrum flow diagram for the
% NRG calculatin. The levels from the same quantum numbers are represendted
% by the same line shape.
% In a plot, the degeneracy is counted at the last iteration to be shown
% (see 'Nshow' option). The degeneracy information is given as '# of
% multiplets (# of states)'. Each legend denotes '[quantum #] (multiplet
% dimension)'.
% The x-axis values for the plot correspond to the NRG iteration for .mat
% files. For example, the energy levels from .._s00.mat (which involves
% only the "impurity") are assigned to Es{1}, and drawn at x = 0.
% 
% 
% < Usage #1 >
% [Etot,Qtot,Qdiff] = plotE (nrgdata [, Options])
%
% Input NRG data file name. 
% 
%   < Input >
%   nrgdata : [char or struct] If char, it is of form 'path/filename' such
%           that the NRG results are saved in 'nrgdata_info.mat' and
%           'nrgdata_##.mat'. If struct, it contains the NRG results.
%
%   < Output >
%   Etot : ['Nshow' x 1 cell] Energy eigenvalues. Each cell (corresponding
%       to each NRG iteration) is another M * 1 cell array whose individual
%       cell corresponds to different quantum # sectors and containes the
%       energy eigenvalues in the iteration.
%   Qtot : ['Nshow' x 1 cell] Quantum #s. Each cell (corresponding to each
%       NRG iteration) is a 2 * 1 cell array where 1st cell contains the
%       quantum # of each symmetry sector, and 2nd cell contains the 
%       multiplet dimension (degeneracy due to non-Abelian symmetry).
%   Qdiff : [numeric] The shifts in quantum numbers per iteration. Each row
%       (say, index n) indicates the shift for the iteration n-1
%       corresponding to Etot{n}. To correctly group the lines in terms of
%       quantum numbers, such shift is compensated, i.e., sum(Qdiff(1:n,:))
%       is subtracted from the quantum number given in Qtot{n}.
%       Note that the first row Qdiff(1,:) for the iteration 0 (which
%       involves the impurity only) is always zero; that is, the quantum
%       number shifts are considered only for the bath sites.
%
%
% < Usage #2 >
%
% [Etot,Qtot,Qdiff] = plotE (Etot,Qtot[, Options])
%
% Input processed data from previous running of plotE.
% NOTE. Usage #2 does not distinguish whether the input is energy or
% entanglement spectrum.
%
%   < Input >
%   Etot, Qtot : Same meaning as the output for usage #1. (But can be further
%          truncated by the option 'Ev'.)
%
%   < Output >
%   Etot, Qtot, Qdiff : Same as the output for usage #1.
% 
%
% < Option (common for Usages #1 and #2) >
% 'beta', .. : [numeric scalar] If it is given, this code draws the flow of
%       the entanglement spectrum, instead of the conventional energy flow
%       diagram. The entanglement spectrum is given by the eigenvalues of the
%       reduced density matrix (= -log(rho)) whose original state is the
%       single-shell construction at the last iteration with fictious inverse
%       temperature 'beta' given by this option.
%       If 'beta' is given in Usage #2, y-axis labels and a title changes
%       to describe entanglement spectrum. In this case, 'Nget' also can be
%       used to appear in title.
%       (Default: [], so draws the conventional energy flow)
% 'dE', .. : [real scalar] Tolerance to determine the emergent degeneracy
%           (not from the non-Abelian symmetry implemented) of eigenstates
%           (Default: 1e-3)
% 'Emax', .. : [real scalar] Upper bound of energy to show (Default: 3)
% 'isFAdded', .. : [logical array] The output of NRG_SL, which is
%       contained in (nrgdata)_info.mat or in the output struct. (For
%       detail, see the description of 'isFAdded' in NRG_SL.m)
%       When the sub-channels in the Wilson chain have different lengths,
%       this optional input is necessary to plot the energy flow diagram.
%       When a whole sub-site of a chain site (associated with a certain
%       sub-channel in the interleaved NRG setting) is not included during
%       the iterative diagonalization, such iteration is skipped in NRG_SL.
%       In this function plotE, the energy levels of such iteration is
%       given by those of the last iteration (with sub-site included)
%       before this skipped iteration, as a workaround. So the result Etot
%       and Qtot change accordingly. This treatment is just ad hoc, but
%       helps the visibility of the flow diagram.
%       Also this optional input is used to determine the proper shift of
%       quantum numbers; see the description of 'Qdiff' above and below.
%       (Default: In the usage #1, retrieve 'isFAdded' from nrgdata if
%       exists. In the usage #2, assume that all the sub-sites are added.) 
% 'Qdiff', .. : [numeric matrix] The shift in the quantum numbers per
%       iteration, which follows the same structure as the output Qdiff.
%       Here this optional input Qdiff overrides the output Qdiff, and also
%       the shift given by this option is compensated in the plot.
%       (Default: no overriding. Qdiff will be determined.)
% 'N2', .. : [integer] Number of iNRG sub-channels.
%       (Default: First, try to read 'N2' from 'nrgdata'_info.mat or the
%       'nrgdata' struct. Then, if 'isFAdded' is given, take its size along
%       the second dimension as 'N2'. Otherwise, use N2 = 1 as default.)
% 'Nget', .. : [integer] Maximal iteration # to collect data from .mat
%       files. For the case of entanglement spectrum flow, 'Nget' is the
%       iteration # where a fictious state is constructed.
%       (Default: For Usage read 'N' in 'nrgdata'_info.mat or the 'nrgdata'
%       struct)
% 'Nshow', .. : [integer] Maximal iteration # to show data. The output
%       'Etot' and 'Qtot' is length of 'Nshow', and the degeneracy and
%       quantum # shift are calculated based on the data at iteratoin
%       'Nshow'.
%       This parameter matters in case of entanglement spectrum, since the
%       entanglement spectum fluctuates within a several steps from the
%       starting state (constructed w/ 'beta' below). For example, the
%       entanglement spectrum at the iteration where the starting state (w/
%       very large 'beta') is constructed would be just a single level w/
%       ground-state degeneracy!
%       (Default: For Usage #1 & energy flow, Nget. For Usage #1 &
%       entanglement spectrum, Nget - 7. See 'Nget' below. 
%       For Usage #2, the length of Etot.)
% 'gcf', .. : [cell array] Option to set the Figure properties.
%       For example, to set the size of figure:
%       'gcf',{'Units','centimeters','Position',[2 2 10 10]}
%       Refer to the MATLAB documentation for Figure properties for
%       details.
%       (Default: [], i.e., no option)
% 'title', .. : [char or cell array] Additonal description to put as the
%       figure title. For single line, give 'title',.. as char array. For
%       multi-line titles, give 'title',.. as the cell array whose elements
%       are char arrays.
%       (Default: [], i.e., no title)
% 'legmax', .. : [numeric scalar] Maximum # of rows to be shown in legend.
%       (Default: 10)
% 'FontSize', .. : [numeric scalar] Font size for the subplot axes.
%       (Default: 12)
% 'noshow' : Suppress showing figure. (Default: shows figure)
% 'save', .. : [char] path/filename for save figures. (Default: [] -> do not save)
%
%
% Written by S.Lee (2015)
% Updated by S.Lee (Oct.13,2017): Now it accepts the struct 'nrgdata' as
%       the NRG result.
% Updated by S.Lee (Feb.21,2018): Become compatible with iNRG.
% Updated by S.Lee (Mar.15,2018): New option 'isFAdded' and 'Qdiff'. Also
%       the convention for Qdiff is somehow changed.
% Updated by S.Lee (Nov.27,2018): To load the results from nrgdata_##.mat,
%       'load3' function is used instead of the standard MATLAB 'load'.
% Updated by S.Lee (Jan.20,2019): Minor fix.
% Updated by S.Lee (Feb.17,2019): Removed option 'panel', and added an
%       option 'N2'.
% Updated by S.Lee (Mar.06,2019): For the iNRG cases, the quantum numbers
%       for individual lines were not properly rounded, so small noise
%       (double-precision error) could differentiate the lines which are of
%       the same quantum numbers. Now this problem is fixed.
% Updated by S.Lee (May 02,2019): Use 'parseQSpace' for parsing the inputs
%       supposed to be QSpace; improved compatibility with AW's
%       NRGWilsonQS.


% % % % for debugging
% try
% % % % 

% default parameters
Ev = 3; 
dE = 1e-2;
oshow = 1;
osave = [];
dstr = [];
Nget = [];
Nshow = [];
beta0 = [];
legmax = 10;
gcfopt = [];
fsize = 12;
isFAdded = [];
Qdiff2 = [];
N2 = [];

if isempty(varargin{1})
    error('ERR: The first input is empty.');
elseif iscell(varargin{1})
    intype = 2; % usage #2
    Etot = varargin{1};
    Qtot = varargin{2};
    varargin(1:2) = [];
elseif ischar(varargin{1})
    intype = 1; % usage #1
    nrgdata = varargin{1};
    varargin(1) = [];
elseif isstruct(varargin{1})
    intype = 3; % usage #1, with 'nrgdata' as the result struct
    nrgdata = varargin{1};
    varargin(1) = [];
elseif (nargin > 1) && isempty(varargin{2})
    error('ERR: The second input is empty.');
elseif ~isempty(varargin{1}) && ~iscell(varargin{1}) && ~ischar(varargin{1}) && ~isstruct(varargin{1})
    error('ERR: The first input is neither cell, char array, nor struct.');
else
    error('ERR: Check usage of plotE.');
end

if nargout > 3
    error('ERR: More than 3 outputs are requested?');
end

while ~isempty(varargin)
    switch varargin{1}
        case 'Emax'
            Ev = varargin{2};
            varargin(1:2) = [];
        case 'dE'
            dE = varargin{2};
            varargin(1:2) = [];
        case 'noshow'
            oshow = 0;
            varargin(1) = [];
        case 'save'
            osave = varargin{2};
            varargin(1:2) = [];
        case 'title'
            dstr = varargin{2};
            varargin(1:2) = [];
        case 'Nget'
            Nget = varargin{2};
            varargin(1:2) = [];
        case 'Nshow'
            Nshow = varargin{2};
            varargin(1:2) = [];
        case 'beta'
            beta0 = varargin{2};
            varargin(1:2) = [];
        case 'legmax'
            legmax = varargin{2};
            varargin(1:2) = [];
        case 'gcf'
            gcfopt = varargin{2};
            varargin(1:2) = [];
        case 'N2'
            N2 = varargin{2};
            varargin(1:2) = [];
        case 'FontSize'
            fsize = varargin{2};
            varargin(1:2) = [];
        case 'isFAdded'
            isFAdded = varargin{2};
            varargin(1:2) = [];
        case 'Qdiff'
            Qdiff2 = varargin{2};
            varargin(1:2) = [];
        otherwise
            disp2(varargin{1});
            error('ERR: Option cannot be interpreted');
    end
end

if any(intype == [1 3]) % usage #1 : read from NRG .mat data
    if isempty(Nget)
        if isempty(beta0) && ~isempty(Nshow)
            Nget = Nshow;
        elseif intype == 3
            Nget = numel(nrgdata.EScale);
            if isempty(isFAdded) && isfield(nrgdata,'isFAdded')
                isFAdded = nrgdata.isFAdded;
            end
        elseif (intype == 1) && exist([nrgdata,'_info.mat'],'file')
            Inrg = load([nrgdata,'_info.mat']);
            Nget = numel(Inrg.EScale);
            if isempty(isFAdded) && isfield(Inrg,'isFAdded')
                isFAdded = Inrg.isFAdded;
            end
        else
            error(['ERR: Usage #1 (read from .mat) failed since ''',nrgdata,'_info.mat'' does not exist.']);
        end
    end
    
    if isempty(Nshow)
        if isempty(beta0)
            Nshow = Nget;
        else
            Nshow = Nget-7;
        end
    end
    
    if isempty(N2)
        if (intype == 1) && isfield(Inrg,'N2')
            N2 = Inrg.N2;
        elseif (intype == 3) && isfield(nrgdata,'N2')
            N2 = nrgdata.N2;
        elseif ~isempty(isFAdded)
            N2 = size(isFAdded,2);
        end
    end        
    
    Etot = cell(Nshow,1); % energy values (minimum = 0, rescaled) 
    Qtot = cell(Nshow,1); % quantum numbers

    % collect energy eigenvalues (& corresponding quantum #s) under threshold
    for itN = (Nget:-1:1)
        if isempty(beta0)
            if intype == 3
                if isempty(nrgdata.HK{itN}) || isempty(nrgdata.HK{itN}.data) % if there is no kept states, use discarded states instead (ex: the last NRG iteration)
                    HK = nrgdata.HT{itN};
                else
                    HK = nrgdata.HK{itN};
                end
            else
                S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'HK');
                HK = S.HK;
                if isempty(HK) || isempty(HK.data) % if there is no kept states, use discarded states instead (ex: the last NRG iteration)
                    S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'HT');
                    if isfield(S,'HT') % for AW's NRGWilsonQS, the convention changed to HD
                        HK = S.HT;
                    else
                        S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'HD');
                        HK = S.HD;
                    end
                end
            end
            HK = parseQSpace(HK);
            
            if itN == 1
                [~,Ieig] = eigQS(HK); % At the 1st iteration ('..._00.mat'), HK is not diagonalized but just a matrix.
                HK = QSpace(Ieig.EK);
            end
        else % entanglement spectrum mode
            if itN == Nget % construct the original state (the mixed state in the 'ground state' subspace)
                if intype == 3
                    if isempty(nrgdata.HK{itN}) || isempty(nrgdata.HK{itN}.data) % if there is no kept states, use discarded states instead (ex: the last NRG iteration)
                        HK = nrgdata.HT{itN};
                        AK = nrgdata.AT{itN};
                    else
                        HK = nrgdata.HK{itN};
                        AK = nrgdata.AK{itN};
                    end
                else
                    S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'HK','AK');
                    HK = S.HK;
                    AK = S.AK;
                    if isempty(HK) || isempty(HK.data) % if there is no kept states, use discarded states instead (ex: the last NRG iteration)
                        S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'HT');
                        if isfield(S,'HT') % for AW's NRGWilsonQS, the convention changed to HD
                            HK = S.HT;
                        else
                            S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'HD');
                            HK = S.HD;
                        end
                        S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'AT');
                        if isfield(S,'AT') % for AW's NRGWilsonQS, the convention changed to HD
                            AK = S.AT;
                        else
                            S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'AD');
                            AK = S.AD;
                        end
                        
                    end
                end
                HK = parseQSpace(HK);
                % Note: HK's name is HK (for code-writing convenience), but the content *will* be 'rho'
                AK = parseQSpace(AK);
                
                zdims = getzdim(HK,1,'-p');
                Ztot = 0;

                for itq = (1:numel(HK.data))
                    HK.data{itq} = exp(-beta0*HK.data{itq});
                    Ztot = Ztot + zdims(itq)*sum(HK.data{itq});
                end

                HK0 = HK/Ztot; % HK0 : (reduced) density matrix. Here it is row vector.
            elseif itN > 1
                if intype == 3
                    AK = nrgdata.AK{itN};
                else
                    S = load3([nrgdata,'_',sprintf('%02i',itN-1),'.mat'],'AK');
                    AK = S.AK;
                end
                if isempty(AK) || isempty(AK.data) % if there is no kept states
                    error(['ERR: No kept states in ''',nrgdata,'_',sprintf('%02i',itN-1),'.mat''?']);
                end
                AK = parseQSpace(AK);
            end
            
            if itN <= Nshow
                % obtain the entanglement spectrum
                if itN == Nget
                    HK = HK0; % row vector; already diagonalized!
                else %if itN > 1
                    [~,Ieig] = eigQS(HK0); % eigenvalues
                    HK = QSpace(Ieig.EK);
                end

                minHK = zeros(numel(HK.data),1);
                for itq = (1:numel(HK.data))
                    HK.data{itq} = -log(HK.data{itq});
                    minHK(itq) = min(HK.data{itq});
                end
                minHK = min(minHK);
                HK = HK-minHK;
            end
            
            % trace out for the next (itN -> itN-1) step
            if itN == Nget
                HK0 = contract(AK,'!1',{AK,'*',diag(HK0)}); 
            elseif itN > 1
                HK0 = contract(AK,'!1',{AK,'*',HK0}); 
            end
        end

        if itN <= Nshow
            EN = HK.data;

            % truncate the data w/ threshold 'Ev'
            for itq = (1:numel(EN))
                EN{itq} = EN{itq}(EN{itq} <= Ev);
                EN{itq} = EN{itq}(:);
            end

            oks = ~cellfun('isempty',EN);
            EN = EN(oks);
            QN = HK.Q{1}(oks,:);
            zdims = getzdim(HK,1,'-p');

            Etot{itN} = EN;
            Qtot{itN} = {QN, zdims(oks)};
        end
    end
else % usage #2 : input energy & quantum # as cell array
    if isempty(Nshow)
        Nshow = numel(Etot);
    end
    
    for itN = (1:Nshow)
        if (iscell(Qtot{itN}) && (size(Qtot{itN}{1},1) ~= numel(Etot{itN}))) || ...
                (~iscell(Qtot{itN}) && (size(Qtot{itN},1) ~= numel(Etot{itN})))
            error(['ERR: # of sectors in energy values and # of quantum numbers do not match: see #',sprintf('%i',itN),'-th element.']);
        end
        
        for itq = (1:numel(Etot{itN}))
            Etot{itN}{itq} = Etot{itN}{itq}(Etot{itN}{itq} <= Ev);
            Etot{itN}{itq} = Etot{itN}{itq}(:);
        end
        oks = ~cellfun('isempty',Etot{itN});
        Etot{itN} = Etot{itN}(oks);
        if iscell(Qtot{itN})
            Qtot{itN} = {Qtot{itN}{1}(oks,:), Qtot{itN}{2}(oks,:)};
        else
            Qtot{itN} = Qtot{itN}(oks,:);
        end
    end
end

if isempty(N2)
    N2 = 1; % default value
end

Npanel = [2 N2]; % Npanel = [# of rows, # of columns]
pNp = prod(Npanel(:)); % period to determine quantum number shifts (ex: for N2 = 1, plotE groups every even iterations and every odd iterations)

if Nshow <= pNp
    error('ERR: # of panels are larger than the maximum iteration to show.');
end

% default setting of isFAdded (if it is not read in the usage #1 or not given in the usage #2)
if isempty(isFAdded)
    isFAdded = true(numel(Etot)-1,1);
end

% find the iterations at which bath sites are not added and assign
% (copy-and-paste) the energy levels of the last iteration (with bath site
% added) before that each skipped iteration, to the skipped iteration.
isNoS = all(~isFAdded,3).';
if any(isNoS(:))
    isNoS = [false;isNoS(:)]; % linearize
    fprintf('WRN: The following iterations are skipped in the iterative diagonalization due to missing/decoupled local space:\n');
    fprintf('%i ',find(isNoS)-1);
    fprintf('\n');

    Etot2 = cell(numel(isNoS),1); Etot2(~isNoS) = Etot;
    Qtot2 = cell(numel(isNoS),1); Qtot2(~isNoS) = Qtot;

    for itN = (2:numel(Etot2))
        if isNoS(itN)
            Etot2(itN) = Etot2(itN-1); Qtot2(itN) = Qtot2(itN-1);
        end
    end

    Etot = Etot2; Qtot = Qtot2;
    fprintf('WRN: Assign the energy levels of the previous iterations to the skipped iterations.\n');
end


if ~isempty(Qdiff2)
    if ~ismatrix(Qdiff2)
        error('ERR: option ''Qdiff'' needs to be 1 or 2 dimensional.');
    elseif iscell(Qtot{1}) && size(Qdiff2,2) ~= size(Qtot{1}{1},2)
        error('ERR: The 2nd dimension (size) of option ''Qdiff'' does not match with quantum number data.');
    elseif ~iscell(Qtot{1}) && size(Qdiff2,2) ~= size(Qtot{1},2)
        error('ERR: The 2nd dimension (size) of option ''Qdiff'' does not match with quantum number data.');
    end
    if size(Qdiff2,1) == 1
        Qdiff2 = repmat(Qdiff2,[numel(Etot),1]);
    end
    if size(Qdiff2,1) ~= numel(Etot)
        error('ERR: The 1st dimension (size) of option ''Qdiff'' does not match with # of iterations in the data.');
    end
    Qdiff = Qdiff2;
else
    FAsum = sum(isFAdded,2); % sum over each supersite of the bath; counts whether flavors are added (1) or not (0) per super-site
    FAsum = permute(bsxfun(@times,FAsum,ones(1,size(isFAdded,2),1)),[2 1 3]);
    FAsum = reshape(FAsum,[size(isFAdded,2)*size(isFAdded,1),size(isFAdded,3)]);
    FAsum((numel(Etot):end),:) = []; % to fit the length of the data, just in case
    
    % last iterations for each set of iterations at which the same combination of flavors are added;
    % +1 due to the iteration 0 for the impurity (which corresponds Etot{1} and Qtot{1})
    isFAsame = [all(FAsum((2:end),:) == FAsum((1:end-1),:),2);false];
    idset = find(~isFAsame) + 1; % necessarily contains numel(Etot) (= the last iteration)

    % find the quantum number shifts
    Qdtot = cell(numel(idset),1);
    for it1 = (1:numel(idset))
        Qdtmp = cell(pNp,1);
        for it2 = (1:pNp)
            if (idset(it1)-it2+1-pNp) > 0 %idset(it1-1)
                Q1 = plotE_getGS(Etot{idset(it1)-it2+1-pNp},Qtot{idset(it1)-it2+1-pNp},dE);
                Q2 = plotE_getGS(Etot{idset(it1)-it2+1},    Qtot{idset(it1)-it2+1},    dE);
                if size(Q1,1) ~= size(Q2,1)
                    fprintf(['WRN: Poor convergence of energy levels at iteration #',sprintf('%i',idset(it1)-it2),': change tolerance dE.\n']);
                end
                Qdtmp{it2} = (mean(Q2,1) - mean(Q1,1))/pNp;
            end
        end
        Qdtmp = cell2mat(Qdtmp);
%         Qdtmp = round(cell2mat(Qdtmp),12); % round up to suppress noise, for fractiaonal quantum number shifts

        if ~isempty(Qdtmp)
            % find most frequent values of shifts
            Qdtmpuniq = unique(Qdtmp,'rows','stable');
            if size(Qdtmpuniq,1) > 1
                Qcnt = zeros(size(Qdtmpuniq,1),1);
                for it2 = (1:size(Qdtmpuniq,1))
                    Qcnt(it2) = sum(all(bsxfun(@eq,Qdtmp,Qdtmpuniq(it2,:)),2));
                end
                [~,maxid] = max(Qcnt);
                Qdtot{it1} = Qdtmpuniq(maxid,:);
            else
                Qdtot{it1} = Qdtmpuniq;
            end
        end
    end

    oktmp = cellfun('isempty',Qdtot);
    idset(oktmp) = []; Qdtot(oktmp) = [];
    
    if ~isempty(Qdtot) && (idset(1) > 1)
        idset = [1;idset];
        Qdtot = [{zeros(1,size(Qdtot{1},2))};Qdtot]; % no quantum number shift for the impurity iteration
        
        for it1 = (2:numel(idset))
            strtmp = sprintf('%.3g, ',Qdtot{it1});
            fprintf(['Detected quantum number shift [',strtmp(1:end-2), ...
                '] per iteration, for iterations (',sprintf('%i:%i',idset(it1-1),idset(it1)-1),')\n']);
            Qdtot{it1} = bsxfun(@times,Qdtot{it1},ones(idset(it1)-idset(it1-1),1));
        end
        Qdiff = cell2mat(Qdtot);
        
    elseif iscell(Qtot{itN})
        Qdiff = zeros(numel(Etot)-1,size(Qtot{end}{1},2));
    else
        Qdiff = zeros(numel(Etot)-1,size(Qtot{end},2));
    end
end

QtotS = Qtot;
Qshift = zeros(1,size(Qdiff,2));
for itN = (1:Nshow)
    Qshift = Qshift + Qdiff(itN,:);
    if iscell(Qtot{itN})
        QtotS{itN}{1} = QtotS{itN}{1} - bsxfun(@times,ones(size(QtotS{itN}{1},1),1),Qshift);
    else
        QtotS{itN} = QtotS{itN} - bsxfun(@times,ones(size(QtotS{itN},1),1),Qshift);
    end
end
% round up, to correct numerical noises; otherwise the same (fractional)
% quantum numbers are regarded to be different due to different ~1e-10 or
% so.
for itN = (1:Nshow)
    if iscell(Qtot{itN})
        QtotS{itN}{1} = round(QtotS{itN}{1},4);
    else
        QtotS{itN} = round(QtotS{itN},4);
    end
end


varargout = cell(1,nargout);
for it = (1:nargout)
    switch it
        case 1
            varargout{it} = Etot;
        case 2
            varargout{it} = Qtot;
        case 3
            varargout{it} = Qdiff;
    end
end
    

if oshow
    % color set (5 default + 2 new)
    lncl = {[0 .447 .741],[.85 .325 .098],[.773 .565 .061], ...
        [.494 .184 .556],[.466 .674 .188],[.301 .745 .933], ...
        [.635 .078 .184]};
    lnst = {'-','--','-.',':'}; % line style
    % These two options are incommensurate; distinguishable are up to 28 lines

    h = figure;
    
    if ~isempty(gcfopt)
        set(h,gcfopt{:});
    else
        set(h,'Units','centimeters','Position',[4 4 (9.5*Npanel(2)) (Npanel(1)*4.5)+3]); % default size of figure window
    end
    
    % To improve readability of RG flow, sort quantum #.
    % List quantum #s and minimum energy for each quantum # sector.
    % Priority : iteration # (descend) -> energy (ascend) -> quantum # (ascend)

    EminQ = cell(Nshow,1); % minimum energy eigenvalues for each sector and corresponding quantum #s
    for itN = (1:Nshow)
        EminQ{itN} = zeros(numel(Etot{itN}),1);
        for itq = (1:numel(Etot{itN}))
            EminQ{itN}(itq) = min(Etot{itN}{itq});
        end
        
        % sort energy (ascend) -> quantum # (ascend) in each iteration
        if iscell(QtotS{itN})
            EminQ{itN} = sortrows([EminQ{itN}, QtotS{itN}{1}, QtotS{itN}{2}]); % also consider multiplet dimension
        else
            EminQ{itN} = sortrows([EminQ{itN}, QtotS{itN}]);
        end
    end
    
    isAllAbel = true(Nshow,1); % check whether the symmetries are all Abelian; if so, the degeneracies will not be printed
    for itN = (1:Nshow)
        if iscell(QtotS{itN}) && any(QtotS{itN}{2} ~= 1)
            isAllAbel(itN) = false;
        end
    end
    isAllAbel = all(isAllAbel);
    
    
    % handles for (mod N) iterations
    hs = cell(Npanel); % for subplots
    hls = cell(Npanel); % for legends
    hxls = cell(Npanel); % for xlabels
    hyls = cell(Npanel); % for ylabels
    ht = []; % for title (of the top subplot panel)
    
    % degeneracy information for (mod N) iterations
    Edeg = cell(pNp,1); % [energy level, multiplet degeneracy (, states degeneracy)]
    
    for itf = (1:pNp) % iterations modulo pNp
        hs{ceil(itf/Npanel(2)),mod(itf-1,Npanel(2))+1} = subplot(Npanel(1),Npanel(2),itf);
        hold on;
        
        Nids = (itf:pNp:Nshow).';
        Quniq = cell2mat(EminQ(flipud(Nids)));
        Quniq = unique(Quniq(:,(2:end)),'rows','stable'); % 'stable' option : keeping order from large iteration # to small iteration #.
        % Quniq: each row indicates the set of quantum numbers. The order
        % of rows are that the energy of lowest level (in the sector
        % labelled by the quantum number) is increasing, from top to
        % bottom.
        legs = cell(size(Quniq,1),1);
        
        if ~isempty(Nids) && ~isempty(Quniq)
            for itQ = (1:size(Quniq,1))
                Etmp = cell(1,numel(Nids));

                for itN2 = (1:numel(Nids))
                    if iscell(QtotS{Nids(itN2)})
                        qok = all(bsxfun(@eq, cell2mat(QtotS{Nids(itN2)}), Quniq(itQ,:)),2);
    %                     qid = find( all( cell2mat(QtotS{Nids(itN2)}) == repmat(Quniq(itQ,:), [size(QtotS{Nids(itN2)}{1},1) 1]), 2) );
                    else
                        qok = all(bsxfun(@eq, QtotS{Nids(itN2)}, Quniq(itQ,:)),2);
    %                     qid = find( all( QtotS{Nids(itN2)} == repmat(Quniq(itQ,:), [size(QtotS{Nids(itN2)},1) 1]), 2) );
                    end
                    if sum(qok) > 1
                        if iscell(QtotS{Nids(itN2)})
                            strtmp = sprintf('%g, ',Quniq(itQ,1:end-1));
                        else
                            strtmp = sprintf('%g, ',Quniq(itQ,:));
                        end
                        strtmp(end-1:end) = [];
                        error(['ERR: In Qtot{',sprintf('%i',Nids(itN2)), ...
                            '}, there are more than one match to quantum number [',strtmp,'].']);
                    end
                    if any(qok)
                        Etmp{itN2} = Etot{Nids(itN2)}{qok};
                    end
    %                 Etmp{itN2} = cell2mat(Etot{Nids(itN2)}(qid(:)));
                end

                maxEnum = max(cellfun('prodofsize',Etmp));

                if maxEnum > 0
                    Etmp2 = nan(numel(Nids)+1,maxEnum);

                    for itN2 = (1:numel(Nids))
                        Etmp2(itN2,(1:numel(Etmp{itN2}))) = sort(Etmp{itN2},'ascend');
                    end

                    Ntmp = [Nids*ones(1,maxEnum); nan(1,maxEnum)];

                    % plot lines
                    plot(Ntmp(:)-1,Etmp2(:), ...
                        'Color',lncl{mod(itQ-1,length(lncl))+1},'LineStyle',lnst{mod(itQ-1,length(lnst))+1},'LineWidth',1);

                    if iscell(QtotS{Nids(1)}) && isAllAbel % if multiplet dimensions are trivially one
                        legs{itQ} = Quniq(itQ,(1:end-1));
                    else
                        legs{itQ} = Quniq(itQ,:);
                    end            
                end
            end

            legs(cellfun('isempty',legs)) = [];
            legs(legmax+1:end) = [];

            legs2 = cell(numel(legs),numel(legs{1})); % each element legs2{n,m} is the conversion of number legs{n}(m) to char
            for it1 = (1:size(legs2,1))
                for it2 = (1:size(legs2,2))
                    legs2{it1,it2} = sprintf('%.4g',legs{it1}(it2));
                    if (legs2{it1,it2}(1) ~= '-') && ...
                            ( (iscell(QtotS{Nids(1)}) && ~isAllAbel && (it2 < size(legs2,2))) || ...
                            ~iscell(QtotS{Nids(1)}) || isAllAbel )
                        legs2{it1,it2} = [' ',legs2{it1,it2}];
                    end
                end
            end
            for it2 = (1:size(legs2,2))
                ndigitmax = max(cellfun('prodofsize',legs2(:,it2)));
                for it1 = (1:size(legs2,1))
                    if (iscell(QtotS{Nids(1)}) && ~isAllAbel && (it2 == size(legs2,2)))
                        legs2{it1,it2} = [' (',legs2{it1,it2},')'];
                    else
                        legs2{it1,it2} = [legs2{it1,it2},repmat(' ',[1 ndigitmax-numel(legs2{it1,it2})])];
                    end
                    if it2 == 1
                        legs2{it1,it2} = ['[',legs2{it1,it2}];
                    end
                    if (iscell(QtotS{Nids(1)}) && ~isAllAbel && (it2 == (size(legs2,2)-1))) || ...
                            ( (~iscell(QtotS{Nids(1)}) || isAllAbel) && (it2 == size(legs2,2)) ) 
                        legs2{it1,it2} = [legs2{it1,it2},']'];
                    end
                end
            end
            for it1 = (1:numel(legs))
                legs{it1} = cell2mat(legs2(it1,:));
            end
            hls{ceil(itf/Npanel(2)),mod(itf-1,Npanel(2))+1} = legend(legs,'Location','eastoutside','Box','off'); %,'FontSize',fsize

            
            % [energy values, multiplet dimension] at the last even/odd iteration
            Edimtmp = cell(numel(Etot{Nids(end)}),1);
            for itq = (1:numel(Etot{Nids(end)}))
                Edimtmp{itq} = [Etot{Nids(end)}{itq}(:), ones(numel(Etot{Nids(end)}{itq}), 1)];

                if iscell(QtotS{Nids(end)})
                    Edimtmp{itq}(:,2) = QtotS{Nids(end)}{2}(itq)*Edimtmp{itq}(:,2);
                end
            end
            Edimtmp = sortrows(cell2mat(Edimtmp));

            % detect degeneracy
            Edeg{itf} = cell(size(Edimtmp,1),1); % [energy value, multiplet degeneracy, state degeneracy]
            cnt = 1;

            while ~isempty(Edimtmp)
                oks = (Edimtmp(:,1) <= (Edimtmp(1,1)+dE));
                Edeg{itf}{cnt} = [Edimtmp(1,1), sum(oks), sum(Edimtmp(oks,2))];
                Edimtmp(oks,:) = [];
                cnt = cnt+1;
            end

            Edeg{itf}(cellfun('isempty',Edeg{itf})) = [];
            Edeg{itf} = cell2mat(Edeg{itf});

            if isAllAbel % all Abelian symmetries
                Edeg{itf}(:,3) = [];
            end
        end
        
        xlim([0 Nshow+0.5]);
        ylim([-0.1 Ev+0.1]);
        set(hs{ceil(itf/Npanel(2)),mod(itf-1,Npanel(2))+1},'FontSize',fsize,'Box','on','LineWidth',1);
        
        if itf == 1
%             titlestr = {dstr,['quantum num. shifted by ',mat2str(-Qdiff),' every site (both bath and impurity)']};
            if iscell(dstr)
                titlestr = [dstr(:);{[]}];
            else
                titlestr = {dstr,[]};
            end
            if ~isempty(beta0)
                titlestr{2} = [titlestr{2},'\beta = ',sprintf('%g',beta0)];
                if ~isempty(Nget)
                    titlestr{2} = [titlestr{2},' @ N = ',sprintf('%i',Nget)];
                end
%             elseif ~isempty(Nget)
%                 titlestr{2} = [titlestr{2},'N = ',sprintf('%i',Nget)];
            end
            titlestr(cellfun('isempty',titlestr)) = [];
            if ~isempty(titlestr)
                if any(cellfun(@(x) ~isempty(strfind(x,'$')), titlestr))
                    ht = title(titlestr,'Interpreter','latex');
                else
                    ht = title(titlestr);
                end
            end
        end
        
        hxls{ceil(itf/Npanel(2)),mod(itf-1,Npanel(2))+1} = ...
            xlabel(['Iteration n \equiv ',sprintf('%i',itf-1),' mod ',sprintf('%i',pNp)]);
        
        if mod(itf-1,Npanel(2)) == 0
            if isempty(beta0)
                hyls{ceil(itf/Npanel(2)),mod(itf-1,Npanel(2))+1} = ylabel('Rescaled energy');
            else
                hyls{ceil(itf/Npanel(2)),mod(itf-1,Npanel(2))+1} = ylabel('log (RDM eigenvalue)');
            end
        end
        
        grid on;
        hold off;
    end
    
    % change the units to 'normalized', to re-organize the panels more easily
    for itf = (1:pNp)
        set(hs{itf},'Units','normalized');
        if ~isempty(hls{itf})
            set(hls{itf},'Units','normalized');
        end
        if ~isempty(hxls{itf})
            set(hxls{itf},'Units','normalized');
        end
        if ~isempty(hyls{itf})
            set(hyls{itf},'Units','normalized');
        end
    end
    if ~isempty(ht)
        set(ht,'Units','normalized');
    end
    
    % size of panels/labels/legends/title
    hss = cell(Npanel);
    hlss = cell(Npanel);
    hxlss = cell(Npanel);
    hylss = cell(Npanel);
    for itf = (1:pNp)
        hss{itf} = get(hs{itf},'Position');
        if ~isempty(hls{itf})
            hlss{itf} = get(hls{itf},'Position');
        end
        if ~isempty(hxls{itf})
            hxlss{itf} = get(hxls{itf},'Extent');
        end
        if ~isempty(hyls{itf})
            hylss{itf} = get(hyls{itf},'Extent');
        end
    end
    if ~isempty(ht)
        hts = get(ht,'Extent');
    end
    
    xsizes = zeros(Npanel(1),1+(4*Npanel(2)));
    ysizes = zeros(1+(4*Npanel(1)),Npanel(2));
    
    for it1 = (1:Npanel(1))
        for it2 = (1:Npanel(2))
            if it2 == 1
                xsizes(it1,1) = hss{it1,it2}(3)*hylss{it1,it2}(3);
                xsizes(it1,2) = hss{it1,it2}(3)*abs(hylss{it1,it2}(1)+hylss{it1,it2}(3));
                xsizes(it1,(6:4:end)) = xsizes(it1,2)*1.25;
                % Extent of hxls{..} is relative to the size of hs{..}
            end
            if ~isempty(hls{it1,it2})
                xsizes(it1,1+(it2*4)) = hlss{it1,it2}(3);
                xsizes(it1,it2*4) = hlss{it1,it2}(1)-hss{it1,it2}(1)-hss{it1,it2}(3);
            end
            if ~isempty(hxls{it1,it2})
                ysizes(1+(it1*4),it2) = hss{it1,it2}(4)*hxlss{it1,it2}(4);
                ysizes(it1*4,it2) = hss{it1,it2}(4)*abs(hxlss{it1,it2}(2)+hxlss{it1,it2}(4));
                % Extent of hxls{..} is relative to the size of hs{..}
            end
        end
    end
    
    ysizes((2:4:end),:) = ysizes((5:4:end),:)*0.5;
    
    if ~isempty(ht)
        it1 = 1; it2 = 1;
        ysizes(it1,it2) = hts(4)*hss{it1,it2}(4);
        ysizes(2,:) = ysizes(2,:)*0.5; % further decrease
    else
        ysizes(2,:) = 0;
    end
    
    xsizes = max(xsizes,[],1); 
    xsizes = [0.01,xsizes,0.01]; % add margins
    xsizes(xsizes<0) = 0; 
    ysizes = max(ysizes,[],2); 
    ysizes = [0.01;ysizes;0.01]; % add margins
    ysizes(ysizes<0) = 0;
    
    xsizes(4:4:end) = (1-sum(xsizes))/Npanel(2);
    ysizes(4:4:end) = (1-sum(ysizes))/Npanel(1);
    
    for it1 = (1:Npanel(1))
        for it2 = (1:Npanel(2))
            set(hs{it1,it2},'Position',[sum(xsizes(1:(4*it2-1))), sum(ysizes((4*it1+1):end)), xsizes(4*it2), ysizes(4*it1)]);
            if ~isempty(hls{it1,it2})
                if hlss{it1,it2}(4) > sum(ysizes((4*it1):(4*it1+2)))
                    fprintf(['WRN: Legend for panel #(',sprintf('%i,%i',it1,it2),') is too long; use option ''legmax'' to decrease its length.\n']);
                end
                set(hls{it1,it2},'Position',[sum(xsizes(1:(4*it2+1))), sum(ysizes((4*it1):end))-hlss{it1,it2}(4), xsizes(4*it2+2), hlss{it1,it2}(4)]);
            end
        end
    end
    
%     if ~isempty(ht)
%         stmp1 = get(hs{1,1},'Position');
%         stmp2 = get(ht,'Extent');
%         % subplot (axes) coordinate -> figure coordinate
%         txtmp = interp1([0 1],[stmp1(1),stmp1(1)+stmp1(3)],[stmp2(1) stmp2(1)+stmp2(3)],'linear','extrap');
%         txtmp = txtmp - (mean(txtmp)-0.5); % center horizontally to figure
%         % figure coordinate -> subplot (axes) coordinate
%         txtmp = interp1([stmp1(1),stmp1(1)+stmp1(3)],[0 1],txtmp(1),'linear','extrap');
%         tytmp = interp1([stmp1(2),stmp1(2)+stmp1(4)],[0 1],1-sum(ysizes(1:2)),'linear','extrap');
%         set(ht,'Position',[txtmp,tytmp,0]);
%     end
    
    % Put labels for degeneracy
    for itf = (1:pNp)
        if ~isempty(Edeg{itf})
            subplot(hs{ceil(itf/Npanel(2)),mod(itf-1,Npanel(2))+1});
            hold on;

            % degeneracy labels
            Edlabs = cell(size(Edeg{itf},1),1);
            x0 = max((itf:pNp:Nshow))-1; % default position

            % first put them without considering overlap
            for itE = (1:size(Edeg{itf},1))
                if size(Edeg{itf},2) > 2
                    strtmp = sprintf('%i(%i)',Edeg{itf}(itE,(2:3)));
                else
                    strtmp = sprintf('%i',Edeg{itf}(itE,2));
                end
                Edlabs{itE} = text(x0,Edeg{itf}(itE,1),strtmp, ...
                    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
                    'BackgroundColor','w','Margin',1e-2,'FontSize',fsize*0.8);
            end

            % shift overlaping labels
            for itE = (2:numel(Edlabs))
                % check overlap
                xtmp = get(Edlabs{itE},'Extent');
                xtmps = cell2mat(cellfun(@(x) get(x,'Extent'),Edlabs(1:itE-1),'UniformOutput',0));
                oks = (xtmps(:,2)+xtmp(:,4) > xtmp(2));
    %             oks = cellfun(@(x) (x(2)+x(4)) > xtmp(2), xtmps);

                if any(oks)
                    xtmps = xtmps(oks,:);
                    [xtmp1,ids] = sort(xtmps(:,1),'ascend');
                    xtmp2 = xtmp1 + xtmps(ids,3);
    %                 xtmp1 = cellfun(@(x) x(1), xtmps);
    %                 xtmp2 = cellfun(@(x) x(1)+x(3), xtmps);
    %                 [xtmp1,ids] = sort(xtmp1,'ascend');
    %                 xtmp2 = xtmp2(ids);
                    xtmp1 = [xtmp1;x0];
                    xtmp2 = [0;xtmp2];

                    idtmp = find(round((xtmp1-xtmp2)*100) >= round(xtmp(3)*100),1,'last');
                    if isempty(idtmp)
                        set(Edlabs{itE},'Position',[x0 xtmp(2)]);
                    else
                        set(Edlabs{itE},'Position',[xtmp1(idtmp) xtmp(2)]);
                    end
                end
            end

            hold off;
        end
    end
    
    if ~isempty(osave)
        saveas(h,osave,'fig');
    end

    if ~logical(oshow)
        close(h);
    end
end

% % % % for debugging
% catch e
%     disp2(getReport(e))
%     disp2('Let''s Debug!');
%     keyboard;
% end

end

function Qres = plotE_getGS (E,Q,dE)
% obtain the quantum numbers for the ground-state subspace.
% E: [cell array] Each cell contains the energy levels for each sector.
% Q: [numeric array] Each row indicates the quantum number for each sector,
%    or [cell array] Q{1} is a numeric array such that each row indicates
%    the quantum number for each sector.
% dE: [scalar] Tolerance to determine the ground-state subspace.
% Qres: [numeric array] Each row indicates the quantum number for each
%   sector of the ground-state subspace

E0 = min(cell2mat(E(:)));
oks = false(numel(E),1);

for itq = (1:numel(E))
    oks(itq) = any(E{itq}(:) <= E0+dE);
end

if iscell(Q)
    Qres = Q{1}(oks,:);
else
    Qres = Q(oks,:);
end

end