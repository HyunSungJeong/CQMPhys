function SampleSpin_XTRG_server(parfn, varargin)

    try % in case of bug

        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [Nkeep, Nsweep, Delta, ChainLen, NumSamples, SampleTau, JobName] = loadvar(partot, ...
                            {'Nkeep', 'Nsweep', 'Delta', 'ChainLen', 'NumSamples', 'SampleTau', 'JobName'}, ...
                                                {[], [], [], [], [], [], []});
        
        disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
        disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
        disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
        disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
        disp2(['LOG_FILE : ', getenv('LOG_FILE')]);
    
        [Sample, rho] = SampleSpin_XTRG(NumSamples, ChainLen, Delta, SampleTau, 'Nkeep', Nkeep, 'Nsweep', Nsweep);

        STG = ['/data/',getenv('USER'),'/XTRG_SpinSamples/',JobName];

        FileName = ['XTRG_SpinSample_NumSamples=', sprintf('%d',NumSamples), '.txt'];

        writematrix(Sample, [STG, filesep, FileName], 'Delimiter', ' ');

        save([STG,'/rho.mat'], 'rho');
        

    catch Err
        disp2(getReport(Err));
        rethrow(Err);
    end
end