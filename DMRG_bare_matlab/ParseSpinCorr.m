function [SpCorr, SpCorr_even, SpCorr_odd] = ParseSpinCorr(Delta, DistMax)

    %% Set parameters
    ChainLen = 500;
    NumSamples = 10000;
    Nkeep = 300;
    Nsweep = 10;

    %% Load MPS
    path = 'C:\Users\hsjun\OneDrive\Physics\Research\CQMPhys\DMRG_bare_matlab\CorrData';

    FilePath = [path, filesep, 'SpCorr_Delta=', sprintf('%.15g',Delta) ,'_ChainLen=', sprintf('%d',ChainLen), '_NumSamples=', ...
                    sprintf('%d',NumSamples), '_Nkeep=', sprintf('%d',Nkeep), '_Nsweep=', sprintf('%d',Nsweep), '.mat'];
    
    Corr = load(FilePath);
    Corr = Corr.Corr;

    %% Parse correlation data
    Ncen = round(ChainLen/2);
    Dist_even = 2:2:DistMax;
    SpCorr_even = zeros(1, numel(Dist_even));
    for it = 1:numel(Dist_even)
        for itS = Ncen-DistMax:Ncen
            SpCorr_even(it) = SpCorr_even(it) + Corr{Dist_even(it)}(itS);
        end
        SpCorr_even(it) = SpCorr_even(it)/(DistMax+1);
    end
    
    Dist_odd = 2:2:DistMax;
    SpCorr_odd = zeros(1, numel(Dist_odd));
    for it = 1:numel(Dist_odd)
        for itS = Ncen-DistMax:Ncen
            SpCorr_odd(it) = SpCorr_odd(it) + Corr{Dist_odd(it)}(itS);
        end
        SpCorr_odd(it) = SpCorr_odd(it)/(DistMax+1);
    end
    
    Dist = 1:DistMax;
    SpCorr = zeros(1, numel(Dist));
    for it = 1:numel(Dist)
        for itS = Ncen-50:Ncen+50
            SpCorr(it) = SpCorr(it) + Corr{Dist(it)}(itS);
        end
        SpCorr(it) = SpCorr(it)/101;
    end
end