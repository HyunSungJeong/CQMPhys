function varargout = SampleSpin(NumSamples, ChainLen, Delta, varargin)
    % <Description>
    % Sequentially samples local spin from the ground state of the 
    % quantum XXZ Heisenberg chain obtained by DMRG(2-site update)
    %
    % <Input>
    % NumSamples : [numeric] Number of spin configurations to sample
    %
    % ChainLen : [numeric] The length of the XXZ Heisenberg chain
    %
    % Delta : [numeric] Spin-z interaction of the XXZ Heisenberg chain
    % 
    % <Option>
    % 'Nkeep', ... : [numeric] Number of states to keep in DMRG (2-site update)
    %                       (Default: 100)
    %
    % 'Nsweep', ... : [numeric] Number of sweeps (pairs of right --> left & left --> right sweeps)
    %                       (Default: 10)
    % 
    % 'update', ... : [char] '1site' or '2site'. 
    %                       If '1site', 1-site update DMRG is performed 
    %                       If '2site', 2-site update DMRG is performed
    %                           (Default: 2-site update)
    %
    % 'SaveFile' : If used, the output is saved as a .txt file instead of a giving it as output variable matlab array output
    %                   (Default: not used)
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

    %% Parse options

    % default options
    Nkeep = 100;
    Nsweep = 10;
    update = '2site';
    SaveFile = false;
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

            case 'update'
                if ~ismember(varargin{2}, {'1site', '2site'})
                    if ischar(varargin{2})
                        error(['ERR: Unknown update option ''',varargin{2},'''']);
                    else
                        error('ERR: Unknown update option');
                    end
                else
                    updateOpt = varargin{2};
                    varargin(1:2) = [];
                end

            case 'SaveFile'
                SaveFile = true;
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

    maxNumCompThreads(3);

    MPO = getMPO(ChainLen, 'XXZ', Delta);   % construct MPO for XXZ Heisenberg chain

    disp(['2-site DMRG with Nkeep = ', sprintf('%d', Nkeep), ' (Init : IterDiag)']);     % Run DMRG
    [~, ~, MPS, ~] = DMRG_GS(MPO, Nkeep, Nsweep, 'MPSinit', 'IterDiag', 'update', update, '-v');
    
    disp('=======================================================================');

    %% Sample from ground state MPS

    % define local operators
    Sz = [1,0;0,-1]/2;      % spin-z operator
    S_r = [0,1;0,0];        % spin raising operator
    S_l = [0,0;1,0];        % spin lowering operator
    I = eye(2);             % identity operator

    % define local projectors
    Proj_up = [1,0;0,0];    % projector to the Sz = 1/2 subspace
    Proj_down = [0,0;0,1];  % projector to the Sz = -1/2 subspace

    % define array to save samples
    Sample = zeros(NumSamples, 15);

    % convert MPS to right canonical form
    [MPS, ~, ~] = getCanonForm(MPS, ChainLen, 'Nkeep', Nkeep); 

    % Direct sampling
    for itS = 1:NumSamples
        Cleft = 1;
        disptime(['Sampling #', sprintf('%d', itS), '/', sprintf('%d', NumSamples)]);

        for itN = 1:ChainLen

            Norm = contract(Cleft, 2, [1,2], getIdentity(Cleft, 2), 2, [1,2]);
            Cleft = Cleft / Norm;
            
            Sz_exp = updateLeft(Cleft, 2, MPS{itN}, Sz, 2, MPS{itN});
            Sz_exp = contract(Sz_exp, 2, [1,2], getIdentity(Sz_exp, 2), 2, [1,2]);
            Prob_up = Sz_exp + 1/2;
           
            if rand() < Prob_up
                Sample(itS, itN) = 1;
                Cleft = updateLeft(Cleft, 2, MPS{itN}, Proj_up, 2, MPS{itN});
                if display
                    disp(['Site #', sprintf('%d', itN), ' : up, Norm = ', sprintf('%.5g', Norm)]);
                end 
            else
                Sample(itS, itN) = 0;
                Cleft = updateLeft(Cleft, 2, MPS{itN}, Proj_down, 2, MPS{itN});
                if display
                    disp(['Site #', sprintf('%d', itN), ' : down, Norm = ', sprintf('%.5g', Norm)]);
                end
            end
        end % itN
    end % itS

    if SaveFile

        dataFolder = [fileparts(mfilename('fullpath')), filesep, 'SampleData'];
        if ~exist(dataFolder, 'dir')
            mkdir(dataFolder);
        end

        if Delta <= -1  % ferromagnetic
            saveFolder = [dataFolder, filesep, 'FM'];
        elseif Delta >= 1   % antiferromagnetic
            saveFolder = [dataFolder, filesep, 'AFM'];
        else    % Luttinger liquid
            saveFolder = [dataFolder, filesep, 'LL'];
        end

        if ~exist(saveFolder, 'dir')
            mkdir(saveFolder);
        end

        FileName = ['Delta=', sprintf('%.15g',Delta), '_L=', sprintf('%d',ChainLen), '_NumSamples=', sprintf('%d', NumSamples), '_Nkeep=', sprintf('%d',Nkeep), '_Nsweep=', sprintf('%d',Nsweep), '.txt'];
        writematrix(Sample, [saveFolder, filesep, FileName], 'Delimiter', ' ');
        
    else
        varargout{1} = Sample;
    end

end