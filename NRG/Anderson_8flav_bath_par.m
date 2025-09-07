function Anderson_8flav_bath_par(varargin)

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
    h_vmem = 250;           % Memory (in GB) to be occupied in clusters
    PE = 32;                % # of cores to be occupied in clusters
    syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
    Nkeep = zeros(1,num_jobs);      % Number of states to be kept
    Lambda = zeros(1,num_jobs);     % NRG discretization parameter
    Hyb = zeros(1,num_jobs);        % Hybridization strength
    T = zeros(1,num_jobs);          % temperature

    for it = 1:num_jobs
        partot(it).PE = PE;

        fprintf(['Nkeep for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).Nkeep = intmp;
        Nkeep(it) = intmp;

        fprintf(['Lambda for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).Lambda = intmp;
        Lambda(it) = intmp;

        fprintf(['Hybridization strength for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).Hyb = intmp;
        Hyb(it) = intmp;

        fprintf(['Temperature for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).T = intmp;
        T(it) = intmp;
    end


    if all(Nkeep == Nkeep(end)) && num_jobs ~= 1
      Nkeeps = ['[',sprintf('%.15g',Nkeep(end)),']x',sprintf('%d',num_jobs)];
    else
      Nkeeps = cellfun(@(x) [sprintf('%.15g',partot(x).Nkeep),','],num2cell(1:num_jobs),'UniformOutput',false);
      Nkeeps = cell2mat(Nkeeps);
      Nkeeps = ['[',Nkeeps(1:end-1),']'];
    end

    if all(Lambda == Lambda(end)) && num_jobs ~= 1
      Lambdas = ['[',sprintf('%.15g',Lambda(end)),']x',sprintf('%d',num_jobs)];
    else
      Lambdas = cellfun(@(x) [sprintf('%.15g',partot(x).Lambda),','],num2cell(1:num_jobs),'UniformOutput',false);
      Lambdas = cell2mat(Lambdas);
      Lambdas = ['[',Lambdas(1:end-1),']'];
    end

    if all(Hyb == Hyb(end)) && num_jobs ~= 1
      Hybs = ['[',sprintf('%.15g',Hyb(end)),']x',sprintf('%d',num_jobs)];
    else
      Hybs = cellfun(@(x) [sprintf('%.15g',partot(x).Hyb),','],num2cell(1:num_jobs),'UniformOutput',false);
      Hybs = cell2mat(Hybs);
      Hybs = ['[',Hybs(1:end-1),']'];
    end
  
    if all(T == T(end)) && num_jobs ~= 1
      Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
    else
      Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
      Ts = cell2mat(Ts);
      Ts = ['[',Ts(1:end-1),']'];
    end
  
    parfn = ['8flav_bath_par_Nkeep=',Nkeeps,'_Lambda',Lambdas,'_Hyb=',Hybs,'_T=',Ts];
    
    dispstruct(partot);
    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
    disp(['Saved to : ', parfn]);
  end