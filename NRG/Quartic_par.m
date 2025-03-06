function Quartic_par (varargin)

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
    h_vmem = 64;            % Memory (in GB) to be occupied in clusters
    PE = 6;                 % # of cores to be occupied in clusters
    syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
    Nkeep = 3000;
    nz = ones(1,num_jobs);
    K_z = zeros(1,num_jobs);        % orbital pseudospin exchange coupling: z-component
    Q = zeros(1,num_jobs);          % quartic-quartic(ddddcccc) interaction strength
    T = zeros(1,num_jobs);          % temperature
  
    for it = (1:num_jobs)
  
      partot(it).nz = nz(it);
      partot(it).PE = PE;
      partot(it).Nkeep = Nkeep;
  
      fprintf(['K_z for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).K_z = intmp;
      K_z(it) = intmp;
  
      fprintf(['Q for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).Q = intmp;
      Q(it) = intmp;
  
      ValidIdx = false;
      while ~ValidIdx
  
        fprintf(['Temperature for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
  
        if isnumeric(intmp) && intmp > 0
          partot(it).T = intmp;
          T(it) = intmp;
          ValidIdx = true;
  
        else
          fprintf('WRN: Invalid Input!\n');
        end
      end
  
      JobName = ['K_z=',sprintf('%.15g',partot(it).K_z), ...
                    '_Q=',sprintf('%.15g',partot(it).Q), ...
                      '_T=',sprintf('%.15g',partot(it).T)];     
  
      partot(it).JobName = JobName;
  
    end
  
    for it = (1:num_jobs)
      if ~exist(['/data/',getenv('USER'),'/Quartic/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)],'dir')
        mkdir(['/data/',getenv('USER'),'/Quartic/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)]);
      end
    end


    if all(K_z == K_z(end)) && num_jobs ~= 1
      K_zs = ['[',sprintf('%.15g',K_z(end)),']x',sprintf('%d',num_jobs)];
    else
      K_zs = cellfun(@(x) [sprintf('%.15g',partot(x).K_z),','],num2cell(1:num_jobs),'UniformOutput',false);
      K_zs = cell2mat(K_zs);
      K_zs = ['[',K_zs(1:end-1),']'];
    end
  
    if all(Q == Q(end)) && num_jobs ~= 1
      Qs = ['[',sprintf('%.15g',Q(end)),']x',sprintf('%d',num_jobs)];
    else
      Qs = cellfun(@(x) [sprintf('%.15g',partot(x).Q),','],num2cell(1:num_jobs),'UniformOutput',false);
      Qs = cell2mat(Qs);
      Qs = ['[',Qs(1:end-1),']'];
    end
  
    if all(T == T(end)) && num_jobs ~= 1
      Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
    else
      Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
      Ts = cell2mat(Ts);
      Ts = ['[',Ts(1:end-1),']'];
    end
  
    parfn = ['Quartic_par_K_z=',K_zs,'_Q=',Qs,'_T=',Ts];
  
    dispstruct(partot);
    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
    disp(['Saved to : ', parfn]);
  end