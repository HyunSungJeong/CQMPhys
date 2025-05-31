function varargout = SampleSpin_XTRG(NumSamples, ChainLen, Delta, SampleTau,  varargin)
    % <Description>
    % Samples local spin-z from the low-temperature density matrix of the 
    % spin-1/2 XXZ Heisenberg chain obtained by XTRG
    %
    % <Input>
    % NumSamples : [numeric] Number of spin configurations to sample
    %
    % ChainLen : [numeric] The length of the XXZ Heisenberg chain
    %
    % Delta : [numeric] Spin-z interaction of the XXZ Heisenberg chain
    %
    % SampleTau : [numeric] The imaginary time(inverse temperature) to sample spin configuration
    % 
    % <Option>
    % 'Nkeep', ... : [numeric] Number of states to keep in XTRG (2-site update)
    %                   (Default: 100)
    %
    % 'Nsweep', ... : [numeric] Number of sweeps (pairs of right --> left & left --> right sweeps) in XTRG
    %                   (Default: 10)
    %
    % 'tau0', ... : [numeric] The imaginary time to start XTRG (must be much smaller than tau0)
    %                   (Default : 0.01)
    %
    % 'SaveFile' : If used, the output is saved as a .txt file instead of a giving it as output variable matlab array output
    %                   (Default: not used)
    %
    % 'getRho' : If used, the density matrix at imaginary time 'SampleTau' is saved as a .mat file
    %                   (Default: not used)
    %
    % 'SaveRho' : If used, the density matrix at imaginary time 'SampleTau' is saved as a .mat file
    % 
    % '-v' : If used, the sampling results are displayed
    %                   (Default: not used)
    %
    % <Output>
    % Sample : [numeric array] Sampled spin configurations. NumSample x ChainLen array.
    %                       Each row is the sampled spin configuration along the chain.
    %                       1 for spin up and 0 for spin down
    %                       Phase label is 2 for ferromagnetic, 3 for Luttinger liquid, and 5 for antiferromagnetic phase
    %                       If 'SaveFile' option is used, this sample is saved in the '\SampleData' folder as a .txt file

    %% Parse inputs
    
    if ~isnumeric(NumSamples)
        error('ERR: ''NumSamples'' must be a number');
    elseif ~isscalar(NumSamples)
        error('ERR: ''NumSamples'' must be a single number');
    elseif mod(NumSamples, 1) ~= 0 || NumSamples <= 0
        error('ERR: ''NumSamples'' must be a positive integer');
    end

    if ~isnumeric(ChainLen)
        error('ERR: ''ChainLen'' must be a number');
    elseif ~isscalar(ChainLen)
        error('ERR: ''ChainLen'' must be a single number');
    elseif mod(ChainLen, 1) ~= 0 || ChainLen < 2
        error('ERR: ''ChainLen'' must be an integer larger than one');
    end

    if ~isnumeric(Delta)
        error('ERR: ''Delta'' must be a number');
    elseif ~isscalar(Delta)
        error('ERR: ''Delta'' must be a single real number');
    end

    if ~isnumeric(SampleTau) || ~isscalar(SampleTau)
        error('ERR: ''SampleTau'' must be a number');
    elseif SampleTau <= 0
        error('ERR: ''SampleTau'' must be positive');
    end

    %% Parse options

    % default options
    Nkeep = 100;
    Nsweep = 10;
    tau0 = 0.01;
    SaveFile = false;
    getRho = false;
    SaveRho = false;
    display = false;

    while ~isempty(varargin)
        switch varargin{1}
            case 'Nkeep'
                if ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
                    error('ERR: ''Nkeep'' must be a positive integer');
                elseif mod(varargin{2},1) ~= 0 || varargin{2} < 1
                    error('ERR: ''Nkeep'' must be a positive integer');
                else
                    Nkeep = varargin{2};
                    varargin(1:2) = [];
                end

            case 'Nsweep'
                if ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
                    error('ERR: ''Nsweep'' must be a positive integer');
                elseif mod(varargin{2},1) ~= 0 || varargin{2} < 1
                    error('ERR: ''Nsweep'' must be a positive integer');
                else
                    Nsweep = varargin{2};
                    varargin(1:2) = [];
                end

            case 'tau0'
                if ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
                    error('ERR: ''tau0'' must be a positive integer');
                elseif mod(varargin{2},1) ~= 0 || varargin{2} < 1
                    error('ERR: ''tau0'' must be a positive integer');
                else
                    if varargin{2} > 0.1
                        disp2('WRN: tau0 = ', sprintf('%.15g', varargin{2}), ' is rather large, and can result in numerical errors in XTRG calculations');
                    end
                    tau0 = varargin{2};
                    varargin(1:2) = [];
                end

            case 'SaveFile'
                SaveFile = true;
                varargin(1) = [];

            case 'getRho'
                getRho = true;
                varargin(1) = [];

            case 'SaveRho'
                SaveRho = true;
                varargin(1) = [];

            case '-v'
                display = true;
                varargin(1) = [];

            otherwise
                if ischar(varargin{1})
                    error(['ERR: unknown input ''',varargin{1},'''']);
                else
                    error('ERR: unknown input');
                end

        end % switch-case
    end % while

    %% Run DMRG

    maxNumCompThreads(4);

    MPO = getMPO(ChainLen, 'XXZ', Delta);   % construct MPO for spin-1/2 XXZ Heisenberg chain

    % Run XTRG
    disp2(['XTRG with Nkeep = ', sprintf('%d', Nkeep), ', Nsweep = ', sprintf('%d', Nsweep)]);
    [rho, ~, ~] = run_XTRG(MPO, tau0, Nkeep, Nsweep, SampleTau);
    
    disp2('=======================================================================');

    %% Sample from low-temp density matrix

    % define local operators
    Sz = [1,0;0,-1]/2;      % spin-z operator
    S_r = [0,1;0,0];        % spin raising operator
    S_l = [0,0;1,0];        % spin lowering operator
    I = eye(2);             % identity operator

    % define local projectors
    Proj_up = [1,0;0,0];    % projector to the Sz = 1/2 subspace
    Proj_down = [0,0;0,1];  % projector to the Sz = -1/2 subspace

    % define array to save samples
    Sample = zeros(NumSamples, ChainLen);

    Tr_R = cell(ChainLen+1, 1);     % right part of the partial traced density matrix
    Tr_R{end} = 1;
    for itN = ChainLen:-1:1
        Tr_R{itN} = contract(rho{itN}, 4, 4, Tr_R{itN+1}, 2, 1);
        Tr_R{itN} = contract(Tr_R{itN}, 3, [1,2], getIdentity(Tr_R{itN}, 1), 2, [1,2]); 
    end

    % Direct sampling
    for itS = 1:NumSamples

        Tr_L = 1;   % left part of the partial traced density matrix

        disptime(['Sampling #', sprintf('%d', itS), '/', sprintf('%d', NumSamples)]);

        for itN = 1:ChainLen

            Tr = contract(Tr_L, 2, 2, rho{itN}, 4, 3);
            Tr = contract(Tr, 4, [2,3], getIdentity(Tr, 2), 2, [2,1]);
            Tr = contract(Tr, 2, 2, Tr_R{itN+1}, 2, 1);
            Tr_L = Tr_L / Tr;
            
            Sz_exp = contract(Tr_L, 2, 2, rho{itN}, 4, 3);
            Sz_exp = contract(Sz_exp, 4, [2,3], Sz, 2, [2,1]);
            Sz_exp = contract(Sz_exp, 2, 2, Tr_R{itN+1}, 2, 1);
            Prob_up = Sz_exp + 1/2;
           
            if rand() < Prob_up
                Sample(itS, itN) = 1;
                Tr_L = contract(Tr_L, 2, 2, rho{itN}, 4, 3);
                Tr_L = contract(Tr_L, 4, [2,3], Proj_up, 2, [2,1]);
                if display
                    disp2(['Site #', sprintf('%d', itN), ' : up', ', Tr = ', sprintf('%.4g',Tr)]);
                end 
            else
                Sample(itS, itN) = 0;
                Tr_L = contract(Tr_L, 2, 2, rho{itN}, 4, 3);
                Tr_L = contract(Tr_L, 4, [2,3], Proj_down, 2, [2,1]);
                if display
                    disp2(['Site #', sprintf('%d', itN), ' : down', ', Tr = ', sprintf('%.4g',Tr)]);
                end
            end
        end % itN
    end % itS
    

    if SaveFile

        dataFolder = [fileparts(mfilename('fullpath')), filesep, 'SampleData'];

        if ~exist(dataFolder, 'dir')
            mkdir(dataFolder);
        end

        Nstep = round(log2(SampleTau)-log2(tau0));
        SampleTau = tau0*(2^Nstep);

        FileName = ['Delta=', sprintf('%.15g',Delta), '_L=', sprintf('%d',ChainLen), '_SampleTau=', sprintf('%.15g', SampleTau), ...
                        '_NumSamples=', sprintf('%d', NumSamples), '_Nkeep=', sprintf('%d',Nkeep), '_Nsweep=', sprintf('%d',Nsweep), '.txt'];
        writematrix(Sample, [dataFolder, filesep, FileName], 'Delimiter', ' ');
        
    else
        varargout{1} = Sample;
        if getRho
            varargout{2} = rho;
        end
    end

    if SaveRho

       rhoFolder = [fileparts(mfilename('fullpath')), filesep, 'rhoData'];

       if ~exist(rhoFolder, 'dir')
            mkdir(rhoFolder);
       end

       FileName = ['Delta=', sprintf('%.15g',Delta), '_L=', sprintf('%d',ChainLen), '_SampleTau=', ...
                        sprintf('%.15g', SampleTau), '_Nkeep=', sprintf('%d',Nkeep), '_Nsweep=', sprintf('%d',Nsweep), '.mat'];

       save([rhoFolder, filesep, FileName], 'rho');


    end

end