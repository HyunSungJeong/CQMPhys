function TCK_par (varargin)

    while numel(varargin) > 0

        if ischar(varargin{1})
            parfn = varargin{1};
            varargin(1) = [];
        else
            disp(varargin{1}); 
            error('ERR: Unknown input/option.');
        end

    end

    fprintf('Number of Jobs?\n');
    num_jobs = input('>>> ');

    syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
    h_vmem = 500;           % Memory (in GB) to be occupied in clusters
    PE = 48;                % # of cores to be occupied in clusters
    syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
    Nkeep = 3000;
    Lambda = 2.5;           % NRG discretization parameter
    getSusc = true;        % whether to compute dynamic susceptibility or not
    getCorr = true;        % whether to compute correlation functions or not
    nz = ones(1,num_jobs);
    J0 = zeros(1,num_jobs);
    T = zeros(1,num_jobs);
    
    for it = (1:num_jobs)

        partot(it).nz = nz(it);
        partot(it).PE = PE;
        partot(it).Nkeep = Nkeep;
        partot(it).getSusc = getSusc;
        partot(it).getCorr = getCorr;
        partot(it).Lambda = Lambda;

        fprintf(['J0 for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).J0 = intmp;
        J0(it) = intmp;

        ValidIdx = false;
        while ~ValidIdx

        fprintf(['Temperature for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');

        if isnumeric(intmp) && intmp > 1e-45 
            partot(it).T = intmp;
            ValidIdx = true;
            T(it) = intmp;
        else
            fprintf('WRN: Invalid Input!\n');
        end
        end

        JobName = ['J0=',sprintf('%.15g',partot(it).J0), ...
                        '_T=',sprintf('%.15g',partot(it).T)];     

        partot(it).JobName = JobName;

    end

    for it = 1:num_jobs

        if ~exist(['/data/',getenv('USER'),'/TCK/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)],'dir')
            mkdir(['/data/',getenv('USER'),'/TCK/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)]);
        end
    end

    if all(J0 == J0(end))
        J0s = ['[',sprintf('%.15g',J0(end)),']x',sprintf('%d',num_jobs)];
    else
        J0s = cellfun(@(x) [sprintf('%.15g',partot(x).J0),','],num2cell(1:num_jobs),'UniformOutput',false);
        J0s = cell2mat(J0s);
        J0s = ['[',J0s(1:end-1),']'];
    end

    if all(T == T(end))
        Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
    else
        Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
        Ts = cell2mat(Ts);
        Ts = ['[',Ts(1:end-1),']'];
    end

    parfn = ['TCK_par_J0=',J0s,'_T=',Ts,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];

    dispstruct(partot);
    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
    disp(['Saved to : ', parfn]);
end