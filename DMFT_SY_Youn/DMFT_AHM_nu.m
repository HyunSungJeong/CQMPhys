function DMFT_AHM_nu(parfn, varargin)

    try % % % in case of bug
        setenv('RC_STORE', go('rcs'));
        partot = job_func_preamble(parfn, varargin{:});
  
        for p = (1:numel(partot))
            [PE, ndfix, T, Lambda, nz, sigmarat, Nkeep, U, r, delta, initS] = loadvar(partot(p), ...
            {'PE', 'ndfix', 'T', 'Lambda', 'nz', 'sigmarat', 'Nkeep', 'U', 'r', 'delta', 'initS'}, ...
            {[], [], [], [], [], [], [], [], [], [], []});
    
            parpool('IdleTimeout', 3500); % % % % %
    
            %% result file path
            resname = go(['glo/AHM/AHM_']);
    
            if ~exist(fileparts(resname), 'dir')
                mkdir(fileparts(resname));
            end
    
            if ~isempty(getenv('SLURM_ARRAY_JOB_ID'))
                resname = [resname, par2str('T', T, 'U', U, 'r', r, 'Lambda', Lambda, 'nz', nz, 'Nk', Nkeep), '_j', getenv('SLURM_ARRAY_JOB_ID'), 't', getenv('SLURM_ARRAY_TASK_ID'), '.mat'];
            else
                resname = [resname, char(datetime('now', 'Format', 'yy.MM.dd_HH-mm-ss')), '.mat'];
            end
    
            disp2(['Result will be saved to: ', resname]);
            nrgdata{1} = cellfun(@(x) go(['data/NRG1_itz=', sprintf('%i', x)]), num2cell(1:nz), 'UniformOutput', false);
            nrgdata{2} = cellfun(@(x) go(['data/NRG2_itz=', sprintf('%i', x)]), num2cell(1:nz), 'UniformOutput', false);
    
            %% DMFT parameters
            Bro = 8;
            % amixes = [0.8, 0.6, 0.3];
            amix = 0.6;
            Ndmft = 50; % % % % %
            cvgth = 1e-3; % % % % % 1e-3
            abstol = 1e-8; % % % % % 1e-8
            ieta = 1i * 10 ^ -4; % % % % % 1i*1e-4
            % iscvgsub = zeros(2, 1);
            RVres = cell(2, 1);
    
            %% NRG parameters
            Etrunc = [];
            ETRUNC = [];
            % N = 3; % % % % %
            N = max(ceil(-2 * log(T / 500) / log(Lambda)) + 10, 20); % % % % %
            Nfit = round(-2 * log(1e-8) / log(Lambda));
            emin = T * 10 ^ -4;
            emax = 1e3;
            estep = 150; % Acont=estep, Adisc=2*estep % % % % %
            Etot = cell(2, 1);
            Qtot = cell(2, 1);
            Qdiff = cell(2, 1);
            %%minTE
    
            %%system parameters
            a = 1;
            t = 1;
            tp = r * t * (1 + delta);
            tm = r * t * (1 - delta);
    
            a1 = [a, a];
            a2 = [a, -a];
            b1 = [pi / a, pi / a];
            b2 = [pi / a, -pi / a];
    
            kpyp = @(kpx) pi / a - abs(kpx);
            kpym = @(kpx) -pi / a + abs(kpx);
            % Kis = [0, 0];
            % Xis = [0, 0];
            % [kplist, kpweight, Nkpfac] = kppatch_DMFT(kpmesh, kpymax);
    
            %% frequency grid
            ocont = getAcont(0, 0, 0, 0, 'emin', emin, 'emax', emax, 'estep', estep); % log mesh
            docon = zeros(numel(ocont), 1);
            docon(1) = 0.5 * (ocont(2) - ocont(1));
            docon(2:end - 1) = 0.5 * (ocont(3:end) - ocont(1:end - 2));
            docon(end) = 0.5 * (ocont(end) - ocont(end - 1));
    
            %% define SEs, RhoV2s, system Hamiltonian...
            RhoV2out = cell(1, 2); % for sublattice
            RhoV2out{1} = zeros(numel(ocont), 2, Ndmft); % Gamma (= hybridization function) // pi : before Broyden
            RhoV2out{2} = RhoV2out{1}; %omega, spin, Ndmft
            RhoV2in = RhoV2out; % fter Broyden // RhoV2 is always diagonal in DCA
            Gfloc = RhoV2out; % cluster(local) Greens function of an f orbital
            SEs = RhoV2out; % self-energy
            nures = zeros(2, 2, Ndmft); %sublattice, spin
            mures = zeros(1, Ndmft);
            mures(1) = U / 2; %half filling initial condition
    
            %% define QSpace objects
            [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1); % two-band model but 1-channel impurity
            [FF, ZF, EF] = setItag('s00', 'op', FF, ZF, IF.E);
            NF = quadOp(FF, FF, []);
            % HU Hamiltonian
            HU = U / 2 * sum(NF) * (sum(NF) - 1);
            A0 = getIdentity(setItag('L00', getvac(EF)), 2, EF, 2, 'K00*', [1 3 2]); % isometry
    
            for ito = (1:numel(FF))
                FHU(ito) = contract(FF(ito), '!1', HU, [1 3 2]) - contract(HU, '!1', FF(ito)); % commutator [FK,HU], to be used for the self-energy trick
            end
    
            %%TODO given self energy, nu
            if isempty(initS)
                SEs{1}(:, :, 1) = -0.1i * 1; % Hilbert transform
                SEs{2}(:, :, 1) = SEs{1}(:, :, 1);
            else
                SEs{1}(:, :, 1) = interp1(initS.ocont, initS.SE{1}, ocont);
                SEs{2}(:, :, 1) = interp1(initS.ocont, initS.SE{2}, ocont);
            end
    
            %% DMFT loop
            for it = (1:Ndmft)
                num_threads_SL(PE);
                muc = mures(it);
        
                fprintf('DMFT iteration No. %d\n', it);
    
                % disp2(intmethod);
                for itsub = 1:2 % sublattice index
                    Gftmp = zeros(numel(ocont), 2);
        
                    for itsp = 1:size(Gfloc{itsub}, 2)  % spin index
                        tic;
            
                        parfor ito = 1:numel(ocont) % parfor % Fourier trafo. of GF. Integrated half of BZ and multiplied by 2, exploiting symmetry
                            Gftmp(ito, itsp) = 2 * a ^ 2 / pi / pi / 4 * integral2(@(kp1, kp2) Gfk(kp1, kp2, ocont(ito), ieta, muc, diag(SEs{itsub}(ito, :, it)), a, t, tp, tm, itsp), -pi / a, pi / a, kpym, kpyp, Method = 'iterated', AbsTol = abstol);
                        end
            
                        toc;
                    end
        
                    Gfloc{itsub}(:, :, it) = Gftmp;
        
                    for itsp = 1:2
                        nures(itsub, itsp, it) = -1 / pi * sum(docon(:) .* imag(Gfloc{itsub}(:, itsp, it)) ./ (1 + exp(ocont(:) / T)));
                    end
        
                    %% Gloc -> Hyb %TODO
                    RhoV2out{itsub}(:, :, it) = -imag((ocont + zeros(1, 2)) + muc - SEs{itsub}(:, :, it) - 1 ./ Gfloc{itsub}(:, :, it)) / pi;
        
                end %itsub
    
                disp2('A');
                disp2(nures(1, :, it));
        
                disp2('B');
                disp2(nures(2, :, it));
        
                %% TODO convergence check // RhoV2old updateHyb check
                % if it == 1
                %   iscvg = false; %Do DMFT at least 1 time
                %   RVres{1} = RhoV2out{1}(:, :, 1); % Do NRG with first-made hyb ftn
                %   RVres{2} = RhoV2out{2}(:, :, 1);
        
                % else
                RhoV2mergein = cell2mat(RhoV2in);       % numel(ocont)x(2*2)xNdmft
                RhoV2mergeout = cell2mat(RhoV2out);     % numel(ocont)x(2*2)xNdmft
                nuressum = reshape(sum(sum(nures, 2), 1), [1, Ndmft]);
                [iscvg, RVmerge, mutmp] = updateHyb(RhoV2mergein(:, :, (1:it)), RhoV2mergeout(:, :, (1:it)), 3, mures(1:it), nuressum(1:it), ndfix, 'cvgth', [cvgth, cvgth], '-v', 'amix', amix, 'Broyden', Bro);
                % RVmerge: hybridization function to be used as the input of the next iteration, A and B sublattices stored in the same variable.
                % end
        
                save(resname, '-v7.3');
                disptime(['Saved data to: ', resname]);
    
                if ~iscvg && (it < Ndmft)
                    RVres{1} = RVmerge(:, 1:2);     % The hybridization function to be used as the input of the next DMFT iteration, sublattice A
                    RVres{2} = RVmerge(:, 3:4);     % The hybridization function to be used as the input of the next DMFT iteration, sublattice B
                    mures(it + 1) = mutmp;
                    RhoV2in{1}(:, :, it + 1) = RVres{1}; % Store hyb ftn used in NRG
                    RhoV2in{2}(:, :, it + 1) = RVres{2};
                elseif iscvg
                    dispbox('DMFT converged.');
                    %TODO reperiodization for DCA
                    break;
                end
    
                %% TODO NRG
                num_threads_SL(8);
                Adiscs = cell(2, 1);
                Aconts = cell(2, 1);
                Adisc2sum = cell(2, 1);
                %Himp Hamiltonian
                epsd =- muc;
                Hepsd = epsd * sum(NF);
    
                for itsub = 1:2
                    Himp = Hepsd + HU + (it < 3) * U * 0.05 * (-1) ^ itsub * (NF(1) - NF(2)) / 2;
                    H0 = contract(A0, '!2*', {Himp, '!1', A0}) +1e-40 * getIdentity(A0, 2); % add "infinitesimal" term to keep all symmetry sectors
                    [ff, gg, dff, dgg] = doZLD(ocont, RVres{itsub}, Lambda, N, nz, 'Nfit', Nfit);
        
                    Adiscs{itsub} = cell(numel(FF) * 3, nz); % discrete data
                    Aconts{itsub} = cell(1, size(Adiscs{itsub}, 1)); % continuous (i.e., broadened) spectral function
        
                    for itz = (1:nz) % for different z shifts
        
                        NRG_SL(nrgdata{itsub}{itz}, H0, A0, Lambda, ff{itz}, FF, ZF, ...
                            gg{itz}, NF, 'Nkeep', Nkeep, 'Etrunc', Etrunc, 'ETRUNC', ETRUNC, ...
                                'dff', dff{itz}, 'dgg', dgg{itz});
            
                        if itz == 1
                            %RG flow
                            [Etot{itsub}, Qtot{itsub}, Qdiff{itsub}] = plotE(nrgdata{itsub}{1}, 'Emax', 7, 'noshow');
                        end
            
                        getRhoFDM(nrgdata{itsub}{itz}, T, '-v', 'Rdiag', true);
                        %TODO
                        [odisc, Adiscs{itsub}(:, itz), sigmak] = getAdisc(nrgdata{itsub}{itz}, [FF(:); FHU(:); FHU(:)], [FF(:); FF(:); FHU(:)], ZF, 'emin', emin, 'emax', emax, 'estep', 2 * estep);
                    end
        
                    for ita = (1:size(Adiscs{itsub}, 1))
                        [~, Aconts{itsub}{ita}] = getAcont(odisc, mean(cell2mat(reshape(Adiscs{itsub}(ita, :), [1 1 nz])), 3), sigmak, T / 5, ...
                            'alphaz', sigmarat / nz, 'smin', log(Lambda) / 15, 'emin', emin, 'emax', emax, 'estep', estep); %mean(,3) : mean along the dim 3
                    end %1=[Fup,Fup], 2=[Fdown,Fdown], 3=[FHUup, Fup], 4=[FHUdown, Fdown], 5=[FHUup, FHUup], 6=[FHUdown, FHUdown]
        
                    Adisc2sum{itsub} = mean(cellfun(@(x) sum(x(:)), Adiscs{itsub}(numel(FF) + (1:numel(FF)), :)), 2); %2by1 matrix
        
                    for itsp = 1:2
                        [SEs{itsub}(:, itsp, it + 1), ~, ~] = SEtrick(ocont, Aconts{itsub}{itsp}, Aconts{itsub}{numel(FF) + itsp}, Adisc2sum{itsub}(itsp), Aconts{itsub}{2 * numel(FF) + itsp}, 'ozin', ocont, 'RhoV2in', RVres{itsub}(:, itsp), 'epsd', epsd);
                    end
        
                    save(resname, '-v7.3');
                    disptime(['Saved data to: ', resname]);
                    disp2('SEtricksave');
                end %itsub
    
            end %DMFT (DCA)
    
        end %partot
  
    catch e
        disp2(getReport(e)); % Report error
    
        if ~(isdeployed || ismcc)
            keyboard
        end
    
        % Save current workspace
        errf = ['Error_', mfilename];
    
        if ~isempty(getenv('SLURM_ARRAY_JOB_ID'))
            errf = [errf, '_j', getenv('SLURM_ARRAY_JOB_ID')];
    
            if ~isempty(getenv('SLURM_ARRAY_TASK_ID'))
            errf = [errf, 't', getenv('SLURM_ARRAY_TASK_ID')];
            end
    
        end
    
        errf = [go('pts'), filesep, errf, char(datetime('now', 'Format', 'yy.MM.dd_HH-mm-ss')), '.mat'];
        save(errf, '-v7.3');
        dispbox('Error! Current workspace is saved in :', errf);
    
        % everything should be before rethrowing error
        rethrow(e);
    end
  
end %function
  
function H = Ham(k, a, t, tp, tm)   % non-interacting Hamiltonian
    H = [-2 * tm * cos((k(1) + k(2)) * a) - 2 * tp * cos((k(1) - k(2)) * a), -2 * t * (cos(k(1) * a) + cos(k(2) * a)); ...
            -2 * t * (cos(k(1) * a) + cos(k(2) * a)), -2 * tp * cos((k(1) + k(2)) * a) - 2 * tm * cos((k(1) - k(2)) * a)
        ];
end
  
function g = Gfk(kpx, kpy, om, ieta, muc, SE, a, t, tp, tm, index)      % Green's function
    g = zeros(1, numel(kpx));
  
    for ikp = 1:numel(kpx)
        k = [kpx(ikp), kpy(ikp)];
        gtmp = inv(((om + ieta + muc) * eye(2) - Ham(k, a, t, tp, tm) - SE));
        size(gtmp);
        g(ikp) = gtmp(index, index);
    end
  
end  