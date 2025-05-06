function Anderson_8flav_par(varargin)

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
    PE = 28;                % # of cores to be occupied in clusters
    syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
    Nkeep = 100;
    nz = ones(1,num_jobs);
    Hyb = zeros(1,num_jobs);        % Hybridization strength
    U = zeros(1,num_jobs);          % Hubbard U
    J = zeros(1,num_jobs);          % Inter-valley Hund coupling
    N0 = zeros(1,num_jobs);         % filling offset
    T = zeros(1,num_jobs);          % temperature
  
    for it = (1:num_jobs)
  
      partot(it).nz = nz(it);
      partot(it).PE = PE;
      partot(it).Nkeep = Nkeep;
  
      fprintf(['Hybridization strength for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).Hyb = intmp;
      Hyb(it) = intmp;
  
      fprintf(['Hubbard U for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).U = intmp;
      U(it) = intmp;

      fprintf(['Hund J for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).J = intmp;
      J(it) = intmp;

      fprintf(['N0 for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).N0 = intmp;
      N0(it) = intmp;
  
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
  
      JobName = ['Hyb=',sprintf('%.15g',partot(it).Hyb), ...
                    '_U=',sprintf('%.15g',partot(it).U), ...
                      '_J=',sprintf('%.15g',partot(it).J), ...
                        '_N0=',sprintf('%.15g',partot(it).N0), ...
                            'T=',sprintf('%.15g',partot(it).T)];     
  
      partot(it).JobName = JobName;
  
    end
  
    for it = (1:num_jobs)
      if ~exist(['/data/',getenv('USER'),'/8flav/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)],'dir')
        mkdir(['/data/',getenv('USER'),'/8flav/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)]);
      end
    end


    if all(Hyb == Hyb(end)) && num_jobs ~= 1
      Hybs = ['[',sprintf('%.15g',Hyb(end)),']x',sprintf('%d',num_jobs)];
    else
      Hybs = cellfun(@(x) [sprintf('%.15g',partot(x).Hyb),','],num2cell(1:num_jobs),'UniformOutput',false);
      Hybs = cell2mat(Hybs);
      Hybs = ['[',Hybs(1:end-1),']'];
    end
  
    if all(U == U(end)) && num_jobs ~= 1
      Us = ['[',sprintf('%.15g',U(end)),']x',sprintf('%d',num_jobs)];
    else
      Us = cellfun(@(x) [sprintf('%.15g',partot(x).U),','],num2cell(1:num_jobs),'UniformOutput',false);
      Us = cell2mat(Us);
      Us = ['[',Us(1:end-1),']'];
    end

    if all(J == J(end)) && num_jobs ~= 1
      Js = ['[',sprintf('%.15g',J(end)),']x',sprintf('%d',num_jobs)];
    else
      Js = cellfun(@(x) [sprintf('%.15g',partot(x).J),','],num2cell(1:num_jobs),'UniformOutput',false);
      Js = cell2mat(Js);
      Js = ['[',Js(1:end-1),']'];
    end

    if all(N0 == N0(end)) && num_jobs ~= 1
      N0s = ['[',sprintf('%.15g',N0(end)),']x',sprintf('%d',num_jobs)];
    else
      N0s = cellfun(@(x) [sprintf('%.15g',partot(x).N0),','],num2cell(1:num_jobs),'UniformOutput',false);
      N0s = cell2mat(N0s);
      N0s = ['[',N0s(1:end-1),']'];
    end
  
    if all(T == T(end)) && num_jobs ~= 1
      Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
    else
      Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
      Ts = cell2mat(Ts);
      Ts = ['[',Ts(1:end-1),']'];
    end
  
    parfn = ['8flav_par_Hyb=',Hybs,'_U=',Us,'_J=',Js,'_N0=',N0s,'_T=',Ts];
  
    dispstruct(partot);
    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
    disp(['Saved to : ', parfn]);
  end