function getSpinCorr_server(parfn, varargin)

    try % in case of bug

        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [MaxDist, range, Nkeep, Nsweep, Delta, ChainLen, NumSamples, JobName] = loadvar(partot, ...
                    {'MaxDist', 'range', 'Nkeep', 'Nsweep', 'Delta', 'ChainLen', 'NumSamples', 'JobName'}, ...
                                            {[], [], [], [], [], [], [], []});
        
        disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
        disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
        disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
        disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
        disp2(['LOG_FILE : ', getenv('LOG_FILE')]);
    
        Corr = getSpinCorr(MaxDist, Delta, ChainLen, Nkeep, 'range', range, 'NumSamples', NumSamples, 'Nsweep', Nsweep);

        CorrName = ['SpCorr_', JobName, '.mat'];

        SavePath = ['/data/hyunsung/DMRG_XXZ/CorrData/', CorrName];

        save(SavePath, 'Corr');

    catch Err
        disp2(getReport(Err));
        rethrow(Err);
    end
end