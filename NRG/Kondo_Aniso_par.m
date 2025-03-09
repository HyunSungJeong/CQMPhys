function Kondo_Aniso_par (varargin)

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
    h_vmem = 95;            % Memory (in GB) to be occupied in clusters
    PE = 10;                % # of cores to be occupied in clusters
    syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
    Nkeep = 3000;
    nz = ones(1,num_jobs);
    J_perp = zeros(1,num_jobs);
    J_z = zeros(1,num_jobs);
    T = zeros(1,num_jobs);
    h = zeros(1,num_jobs);
    %h = 0:0.1:1;
    
    for it = (1:num_jobs)
  
      partot(it).nz = nz(it);
      partot(it).PE = PE;
      partot(it).Nkeep = Nkeep;
      partot(it).h = h(it);
      
      fprintf(['J_perp for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).J_perp = intmp;
      J_perp(it) = intmp;
  
      fprintf(['J_z for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).J_z = intmp;
      J_z(it) = intmp;
  
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
  
      if h(it) == 0
        JobName = ['J_perp=',sprintf('%.15g',partot(it).J_perp), ...
                    '_J_z=',sprintf('%.15g',partot(it).J_z), ...
                      '_T=',sprintf('%.15g',partot(it).T)];
      else
        JobName = ['J_perp=',sprintf('%.15g',partot(it).J_perp), ...
                    '_J_z=',sprintf('%.15g',partot(it).J_z), ...
                      '_T=',sprintf('%.15g',partot(it).T), ...
                        '_h=',sprintf('%.15g',partot(it).h)];
      end     
  
      partot(it).JobName = JobName;
  
    end
  
    for it = (1:num_jobs)
      if ~exist(['/data/',getenv('USER'),'/Kondo_Aniso/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)],'dir')
        mkdir(['/data/',getenv('USER'),'/Kondo_Aniso/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)]);
      end
    end

    if all(J_perp == J_perp(end)) && num_jobs ~= 1
      J_perps = ['[',sprintf('%.15g',J_perp(end)),']x',sprintf('%d',num_jobs)];
    else
      J_perps = cellfun(@(x) [sprintf('%.15g',partot(x).J_perp),','],num2cell(1:num_jobs),'UniformOutput',false);
      J_perps = cell2mat(J_perps);
      J_perps = ['[',J_perps(1:end-1),']'];
    end
  
    if all(J_z == J_z(end)) && num_jobs ~= 1
      J_zs = ['[',sprintf('%.15g',J_z(end)),']x',sprintf('%d',num_jobs)];
    else
      J_zs = cellfun(@(x) [sprintf('%.15g',partot(x).J_z),','],num2cell(1:num_jobs),'UniformOutput',false);
      J_zs = cell2mat(J_zs);
      J_zs = ['[',J_zs(1:end-1),']'];
    end
  
    if all(T == T(end)) && num_jobs ~= 1
      Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
    else
      Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
      Ts = cell2mat(Ts);
      Ts = ['[',Ts(1:end-1),']'];
    end

    if all(h == h(end)) && num_jobs ~= 1
      hs = ['[',sprintf('%.15g',h(end)),']x',sprintf('%d',num_jobs)];
    else
      hs = cellfun(@(x) [sprintf('%.15g',partot(x).h),','],num2cell(1:num_jobs),'UniformOutput',false);
      hs = cell2mat(hs);
      hs = ['[',hs(1:end-1),']'];
    end
  
    if all(h == 0)
      parfn = ['Kondo_Aniso_par_J_perp=',J_perps,'_J_z=',J_zs,'_T=',Ts];
    else
      parfn = ['Kondo_Aniso_par_J_perp=',J_perps,'_J_z=',J_zs,'_T=',Ts,'_h=',hs];
    end
  
    dispstruct(partot);
    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
    disp(['Saved to : ', parfn]);
  end