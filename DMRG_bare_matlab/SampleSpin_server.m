function SampleSpin_server(parfn, varargin)

    try % in case of bug

        setenv('RC_STORE', go('rcs'));      % '/project/hyunsung/RCStore'

        partot = job_func_preamble(parfn, varargin{:});
    
        [Nkeep, Nsweep, Delta, ChainLen, NumSamples, JobName] = loadvar(partot, ...
                            {'Nkeep', 'Nsweep', 'Delta', 'ChainLen', 'NumSamples', 'JobName'}, ...
                                                {[], [], [], [], [], []});
        
        disp2(['SLURM_JOB_NAME : ', getenv('SLURM_JOB_NAME')]);
        disp2(['SLURM_JOB_ID : ', getenv('SLURM_JOB_ID')]);
        disp2(['SLURM_ARRAY_JOB_ID : ', getenv('SLURM_ARRAY_JOB_ID')]);
        disp2(['SLURM_ARRAY_TASK_ID : ', getenv('SLURM_ARRAY_TASK_ID')]);
        disp2(['LOG_FILE : ', getenv('LOG_FILE')]);
    
        [Sample, MPS_GS] = SampleSpin(NumSamples, ChainLen, Delta, 'Nkeep', Nkeep, 'Nsweep', Nsweep, 'getMPS');

        STG = ['/data/',getenv('USER'),'/DMRG_SpinSamples/',JobName, '_Nkeep=', sprintf('%.15g',Nkeep), '_Nsweep=', sprintf('%.15g',Nsweep)];

        FileName = 'DMRG_SpinSample.txt';

        writematrix(Sample, [STG, filesep, FileName], 'Delimiter', ' ');

        save([STG,'/MPS_GS.mat'], 'MPS_GS');
        

    catch Err
        disp2(getReport(Err));
        rethrow(Err);
    end
end