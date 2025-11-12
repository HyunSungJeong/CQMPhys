function ToAH_par(varargin)

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
    h_vmem = 500;           % Memory (in GB) to be occupied in clusters
    PE = 48;                % # of cores to be occupied in clusters
    syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
    Nkeep = 3000;
    Lambda = 4;             % NRG discretization parameter 
    N_SuscIter = 2;         % number of iterations to compute susceptibilities
    getSusc = true;         % whether to compute dynamic susceptibilities or not
    nz = ones(1,num_jobs);          % number of z-shifts
    Hyb = zeros(1,num_jobs);        % hybridization strength
    U = zeros(1,num_jobs);          % Hubbard U
    J = zeros(1,num_jobs);          % spin Hund coupling
    J_L = zeros(1,num_jobs);        % orbital Hund coupling
    mu = zeros(1,num_jobs);         % chemical potential
    T = zeros(1,num_jobs);          % temperature

    option = 0;
    while ~ismember(option, [1,2])
      fprintf('1 for mu sweep, 2 for parameter sweep\n');
      option = input('>>> ');
      if ~ismember(option, [1,2])
        fprintf('WRN: Invalid Input!\n');
      end
    end 

    if option == 1  
      
      fprintf('Hybridization strength for all jobs\n');
      Hyb = input('>>> ');
    
      fprintf('Hubbard U for all jobs\n');
      U = input('>>> ');

      fprintf('Spin Hund coupling J for all jobs\n');
      J = input('>>> ');

      fprintf('Orbital Hund coupling J_L for all jobs\n');
      J_L = input('>>> ');
    
      ValidIdx = false;
      while ~ValidIdx
    
        fprintf('Temperature for all jobs\n');
        intmp = input('>>> ');
    
        if isnumeric(intmp) && intmp > 0
          T = intmp;
          ValidIdx = true;
    
        else
          fprintf('WRN: Invalid Input!\n');
        end
      end

      for it = 1:num_jobs

        partot(it).nz = nz(it);
        partot(it).PE = PE;
        partot(it).Nkeep = Nkeep;
        partot(it).Lambda = Lambda;
        partot(it).getSusc = getSusc;
        partot(it).Hyb = Hyb;
        partot(it).U = U;
        partot(it).J = J;
        partot(it).J_L = J_L;
        partot(it).T = T;
        partot(it).N_SuscIter = N_SuscIter;

        fprintf(['Chemical potential for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).mu = intmp;
        mu(it) = intmp;

        JobName = ['Hyb=',sprintf('%.15g',partot(it).Hyb), ...
                    '_U=',sprintf('%.15g',partot(it).U), ...
                      '_J=',sprintf('%.15g',partot(it).J), ...
                        '_J_L=',sprintf('%.15g',partot(it).J_L), ...
                          '_T=',sprintf('%.15g',partot(it).T), ...
                            '_mu=',sprintf('%.15g',partot(it).mu)];     
    
        partot(it).JobName = JobName;

      end

    else    % parameter sweep

      fprintf('Chemical potential for all jobs\n');
      intmp = input('>>> ');
      mu = intmp*ones(1,num_jobs);
  
      for it = 1:num_jobs
    
        partot(it).nz = nz(it);
        partot(it).PE = PE;
        partot(it).Nkeep = Nkeep;
        partot(it).Lambda = Lambda;
        partot(it).getSusc = getSusc;
        partot(it).mu = mu(it);
        partot(it).N_SuscIter = N_SuscIter;
    
        fprintf(['Hybridization strength for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).Hyb = intmp;
        Hyb(it) = intmp;
    
        fprintf(['Hubbard U for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).U = intmp;
        U(it) = intmp;

        fprintf(['Spin Hund coupling J for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).J = intmp;
        J(it) = intmp;

        fprintf(['Orbital Hund coupling J_L for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
        partot(it).J_L = intmp;
        J_L(it) = intmp;
    
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
                          '_J_L=',sprintf('%.15g',partot(it).J_L), ...
                              '_T=',sprintf('%.15g',partot(it).T), ...
                                '_mu=',sprintf('%.15g',partot(it).mu)];     
    
        partot(it).JobName = JobName;
    
      end
    
    end   % option
  
    for it = 1:num_jobs
      if ~exist(['/data/',getenv('USER'),'/ToAH/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)],'dir')
        mkdir(['/data/',getenv('USER'),'/ToAH/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)]);
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

    if all(J_L == J_L(end)) && num_jobs ~= 1
      J_Ls = ['[',sprintf('%.15g',J_L(end)),']x',sprintf('%d',num_jobs)];
    else
      J_Ls = cellfun(@(x) [sprintf('%.15g',partot(x).J_L),','],num2cell(1:num_jobs),'UniformOutput',false);
      J_Ls = cell2mat(J_Ls);
      J_Ls = ['[',J_Ls(1:end-1),']'];
    end
  
    if all(T == T(end)) && num_jobs ~= 1
      Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
    else
      Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
      Ts = cell2mat(Ts);
      Ts = ['[',Ts(1:end-1),']'];
    end

    if all(mu == mu(end)) && num_jobs ~= 1
      mus = ['[',sprintf('%.15g',mu(end)),']x',sprintf('%d',num_jobs)];
    else
      mus = cellfun(@(x) [sprintf('%.15g',partot(x).mu),','],num2cell(1:num_jobs),'UniformOutput',false);
      mus = cell2mat(mus);
      mus = ['[',mus(1:end-1),']'];
    end
  
    if getSusc
      parfn = ['Susc_ToAH_par_Hyb=',Hybs,'_U=',Us,'_J=',Js,'_J_L=',J_Ls,'_T=',Ts,'_mu=',mus,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];
    else
      if option == 1
        parfn = ['MuSweep_ToAH_par_Hyb=',Hybs,'_U=',Us,'_J=',Js,'_J_L=',J_Ls,'_T=',Ts,'_mu=',mus,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];
      else
        parfn = ['Entropy_ToAH_par_Hyb=',Hybs,'_U=',Us,'_J=',Js,'_J_L=',J_Ls,'_T=',Ts,'_mu=',mus,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];
      end
    end
  
    dispstruct(partot);
    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
    disp(['Saved to : ', parfn]);
  end