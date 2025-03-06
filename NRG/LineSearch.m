function LineSearch(parfn,varargin)

    setenv('RC_STORE', go('rcs'));      % '/home/hyunsung/RCStore'

    partot = job_func_preamble(parfn, varargin{:});

    [PE, Bi_Par, T, J0, K0, I0, hi, lo, JobName] = loadvar(partot, ...
        {'PE', 'Bi_Par', 'T', 'J0', 'K0', 'I0', 'hi', 'lo', 'JobName'}, ...
            {[], [], [], [], [], [], [], [], []});

    strtmp = cell(5,1);
    strtmp{1} = ['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')];
    strtmp{2} = ['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')];
    strtmp{3} = ['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')];
    strtmp{4} = ['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')];
    strtmp{5} = ['LOG_FILE : ', getenv('LOG_FILE')];
    dispbox('-width',130,strtmp{:});

    Storage_dir = ['/data/',getenv('USER'),'/LineSearch/',Bi_Par,'_search/',JobName];

    MaxSep = 1e-2;      % Maximum separation of energy levels that are regarded as degenerate
    Thres = 1e-13;       % Maximum difference of system parameters that are regarded as different
    NumComp = 3;       % Number of lowest degenerate multiplets to be compared
    num_threads_SL(8);

    disp2('');
    disp2('*Caution: This algorithm is guaranteed to work');
    disp2('only if there are no two disconnected intervals with the same phase');
    disp2('');
    disp2('Otherwise, this algorithm must be repeatedly applied');
    disp2('to intervals that are suspected to contain phase boundaries');
    disp2('');
    disptime('Start Search');
    strtmp = cell(3,1);
    strtmp{1} = ['Searching Interval ',Bi_Par,' = [',sprintf('%.15g',lo),',',sprintf('%.15g',hi),']'];

    if Bi_Par == 'J0'
        strtmp{2} = ['K0 = ',sprintf('%.15g',K0),', I0 = ',sprintf('%.15g',I0)];
    elseif Bi_Par == 'K0'
        strtmp{2} = ['J0 = ',sprintf('%.15g',J0),', I0 = ',sprintf('%.15g',I0)];
    elseif Bi_Par == 'I0'
        strtmp{2} = ['J0 = ',sprintf('%.15g',J0),', K0 = ',sprintf('%.15g',K0)];
    end

    strtmp{3} = ['MaxSep = ',sprintf('%.15g',MaxSep),', Thres = ',sprintf('%.15g',Thres)];
    dispbox('-width',75,strtmp{:});
    disp2('');

    % NRG parameters
    % ***************************************************
    %% Lambda and Nkeep must be syncronized with BiSearch
    % ***************************************************
    Lambda = 4;
    N = max(ceil(-2*log(T/100)/log(Lambda))+8,20);
    Nkeep = 5000;        

    [ff, gg] = doZLD([-1;1],[1;1],Lambda,N,1,'Nfit',round(-2*log(1e-8)/log(Lambda)));

    if Bi_Par == 'J0'
        [DegMul_lo, Etot_lo, Qtot_lo] = DegMul(lo,K0,I0,ff,Nkeep,Lambda,MaxSep);
        DegMul_lo = DegMul_lo(1:NumComp);
        [DegMul_hi, Etot_hi, Qtot_hi] = DegMul(hi,K0,I0,ff,Nkeep,Lambda,MaxSep);
        DegMul_hi = DegMul_hi(1:NumComp);
    elseif Bi_Par == 'K0'
        [DegMul_lo, Etot_lo, Qtot_lo] = DegMul(J0,lo,I0,ff,Nkeep,Lambda,MaxSep);
        DegMul_lo = DegMul_lo(1:NumComp);
        [DegMul_hi, Etot_hi, Qtot_hi] = DegMul(J0,hi,I0,ff,Nkeep,Lambda,MaxSep);
        DegMul_hi = DegMul_hi(1:NumComp);
    elseif Bi_Par == 'I0'
        [DegMul_lo, Etot_lo, Qtot_lo] = DegMul(J0,K0,lo,ff,Nkeep,Lambda,MaxSep);
        DegMul_lo = DegMul_lo(1:NumComp);
        [DegMul_hi, Etot_hi, Qtot_hi] = DegMul(J0,K0,hi,ff,Nkeep,Lambda,MaxSep);
        DegMul_hi = DegMul_hi(1:NumComp);
    end

    save([Storage_dir,'/Etot_',Bi_Par,'=',sprintf('%.15g',lo),'.mat'],'Etot_lo');
    save([Storage_dir,'/Qtot_',Bi_Par,'=',sprintf('%.15g',lo),'.mat'],'Qtot_lo');
    save([Storage_dir,'/Etot_',Bi_Par,'=',sprintf('%.15g',hi),'.mat'],'Etot_hi');
    save([Storage_dir,'/Qtot_',Bi_Par,'=',sprintf('%.15g',hi),'.mat'],'Qtot_hi');

    if isequal(DegMul_hi,DegMul_lo)

        strtmp = cell(6,1);
        strtmp{1} = ['Phase at ',Bi_Par,' = ',sprintf('%.15g',lo),' and ',Bi_Par,' = ',sprintf('%.15g',hi),' is the same'];
        strtmp{2} = 'Number of degenerate multiplets at both points :';
        tmp = cellfun(@(x) [sprintf('%.15g',DegMul_lo(x)),' '],num2cell(1:numel(DegMul_lo)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{3} = ['[ ',tmp,']'];
        strtmp{4} = ' ';
        strtmp{5} = 'Probing phase at';
        tmp = lo + (hi-lo)*(1:9)/10;
        tmp = cellfun(@(x) [sprintf('%.15g',tmp(x)),' '],num2cell(1:numel(tmp)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{6} = [Bi_Par,' = [ ',tmp,']'];
        dispbox('-width',130,strtmp{:});

        DegMuls = cell(9,1);

        for it = (1:9)

            par = lo + (hi-lo)*it/10;
            if Bi_Par == 'J0'
                [DegMul_par, Etot, Qtot] = DegMul(par,K0,I0,ff,Nkeep,Lambda,MaxSep);
                DegMul_par = DegMul_par(1:NumComp);
            elseif Bi_Par == 'K0'
                [DegMul_par, Etot, Qtot] = DegMul(J0,par,I0,ff,Nkeep,Lambda,MaxSep);
                DegMul_par = DegMul_par(1:NumComp);
            elseif Bi_Par == 'I0'
                [DegMul_par, Etot, Qtot] = DegMul(J0,K0,par,ff,Nkeep,Lambda,MaxSep);
                DegMul_par = DegMul_par(1:NumComp);
            end

            strtmp = cell(2,1);
            strtmp{1} = ['Numbers of degenerate multiplets at ',Bi_Par,' = ',sprintf('%.15g',par),' :'];
            tmp = cellfun(@(x) [sprintf('%.15g',DegMul_par(x)),' '],num2cell(1:numel(DegMul_par)),'UniformOutput',false);
            tmp = cell2mat(tmp);
            strtmp{2} = ['[ ',tmp,']'];
            dispbox('-width',130,strtmp{:});

            DegMuls{it} = DegMul_par;

            save([Storage_dir,'/Etot_',Bi_Par,'=',sprintf('%.15g',par),'.mat'],'Etot');
            save([Storage_dir,'/Qtot_',Bi_Par,'=',sprintf('%.15g',par),'.mat'],'Qtot');
            save([Storage_dir,'/Num_DegMul_',Bi_Par,'=',sprintf('%.15g',par),'.mat'],'DegMul_par');

        end

        strtmp = cell(11,1);
        strtmp{1} = 'Numbers of degenerate multiplets at ...';
        strtmp{2} = ' ';
        for it = (1:9)

            DegMul_par = DegMuls{it};
            tmp = cellfun(@(x) [sprintf('%.15g',DegMul_par(x)),' '],num2cell(1:numel(DegMul_par)),'UniformOutput',false);
            tmp = cell2mat(tmp);
            strtmp{it+2} = [Bi_Par,' = ',sprintf('%.15g',lo + (hi-lo)*it/10),' : [ ',tmp,']'];
        end
        dispbox('-width',130,strtmp{:});

    else

        [Boundaries, Num_DegMuls] = DFS_Tree(Bi_Par,T,J0,K0,I0,lo,hi,DegMul_lo,DegMul_hi,MaxSep,Thres,Storage_dir,'NumComp',NumComp);
        
        Boundaries = cat(2,lo,Boundaries,hi);
        disp2('');
        disptime('Search Complete');
        strtmp = cell(4+numel(Num_DegMuls),1);
        strtmp{1} = 'Boundaries :';
        tmp = cellfun(@(x) [sprintf('%.15g',Boundaries(x)),' '],num2cell(1:numel(Boundaries)),'UniformOutput',false);
        tmp = cell2mat(tmp);
        strtmp{2} = ['[ ',tmp,']'];
        strtmp{3} = ' ';
        strtmp{4} = 'Numbers of degenerate multiplets at each interval :';
        for it = (1:numel(Num_DegMuls))
            tmp1 = Num_DegMuls{it};
            tmp2 = cellfun(@(x) [sprintf('%.15g',tmp1(x)),' '],num2cell(1:numel(tmp1)),'UniformOutput',false);
            tmp2 = cell2mat(tmp2);
            strtmp{4+it} = ['[ ',tmp2,']'];
        end
        dispbox('-width',130,strtmp{:});

        save([storage_dir,'/Boundaries.mat'],'Boundaries');
        save([storage_dir,'/Num_DegMuls.mat'],'Num_DegMuls');
    end

    LOG_FilePath = getenv('LOG_FILE');
    ERR_FilePath = [Log_FilePath(1:end-4), '_err', Log_FilePath(end-3:end)];

    if Bi_Par == 'J0'
        movefile(LOG_FilePath, [getenv('STORAGE_DIR'),'/LineSearch/J0_search/',JobName]);
        movefile(ERR_FilePath, [getenv('STORAGE_DIR'),'/LineSearch/J0_search/',JobName]);
    elseif Bi_Par == 'K0'
        movefile(LOG_FilePath, [getenv('STORAGE_DIR'),'/LineSearch/K0_search/',JobName]);
        movefile(ERR_FilePath, [getenv('STORAGE_DIR'),'/LineSearch/K0_search/',JobName]);
    elseif Bi_Par == 'I0'
        movefile(LOG_FilePath, [getenv('STORAGE_DIR'),'/LineSearch/I0_search/',JobName]);
        movefile(ERR_FilePath, [getenv('STORAGE_DIR'),'/LineSearch/I0_search/',JobName]);
    end
    
end