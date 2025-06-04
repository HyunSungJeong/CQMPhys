function getSpinCorr_par(varargin)

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
    h_vmem = 30;            % Memory (in GB) to be occupied in clusters
    PE = 3;                 % # of cores to be occupied in clusters
    Nkeep = 300;
    Nsweep = 10;
    NumSamples = 10000;
    Delta = ones(1,num_jobs);

    fprintf('MaxDist for all jobs?\n');
    MaxDist = input('>>> ');

    fprintf('Range for all jobs (in the form of [min,max])?\n');
    range = input('>>> ');

    fprintf('ChainLen for all jobs?\n');
    ChainLen = input('>>> ');
  
    for it = (1:num_jobs)

      partot(it).PE = PE;
      partot(it).Nkeep = Nkeep;
      partot(it).Nsweep = Nsweep;
      partot(it).ChainLen = ChainLen;
      partot(it).NumSamples = NumSamples;
      partot(it).MaxDist = MaxDist;
      partot(it).range = range;
  
      fprintf(['Delta for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).Delta = intmp;
      Delta(it) = intmp;
 
      JobName = ['Delta=',sprintf('%.15g',partot(it).Delta), ...
                    '_ChainLen=',sprintf('%.15g',ChainLen), ...
                      '_NumSamples=',sprintf('%.15g',NumSamples), ...
                        '_Nkeep=',sprintf('%.15g',Nkeep), ...
                          '_Nsweep=',sprintf('%.15g',Nsweep)];     
  
      partot(it).JobName = JobName;
  
    end

    if all(Delta == Delta(end)) && num_jobs ~= 1
      Deltas = ['[',sprintf('%.15g',Delta(end)),']x',sprintf('%d',num_jobs)];
    else
      Deltas = cellfun(@(x) [sprintf('%.15g',partot(x).Delta),','],num2cell(1:num_jobs),'UniformOutput',false);
      Deltas = cell2mat(Deltas);
      Deltas = ['[',Deltas(1:end-1),']'];
    end
  
    parfn = ['DMRG_SpinCorr_par_Delta=', Deltas, '_NumSamples=', sprintf('%d', NumSamples), '_ChainLen=', sprintf('%d', ChainLen), '_Nkeep=', sprintf('%d', Nkeep)];
  
    dispstruct(partot);
    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
    disp(['Saved to : ', parfn]);
  end