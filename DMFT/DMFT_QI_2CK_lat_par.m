function DMFT_QI_2CK_lat_par (varargin)

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
    Nkeep = 3000;
    Lambda = 4;
    sigmarat = 1.5;
    nz = 1;
    emax = 1e3;         % Maximum absolute value of frequency grid
    estep = 250;        % # of steps to increase frequency per decade(x10), used in getAdisc and getAcont
    ndfix = zeros(1, num_jobs);       % target filling
    U = zeros(1, num_jobs);           % Coulomb repulsion energy
    V = zeros(1, num_jobs);           % inter-cell hopping amplitude V
    t_0 = zeros(1, num_jobs);         % intra-cell hopping amplitude t_0
    phi_div_pi = zeros(1, num_jobs);  % distortion angle phi
    T = zeros(1, num_jobs);           % temperature
  
    for it = (1:num_jobs)
  
      partot(it).PE = PE;
      partot(it).Nkeep = Nkeep;
      partot(it).Lambda = Lambda;
      partot(it).sigmarat = sigmarat;
      partot(it).nz = nz;
      partot(it).emax = emax;
      partot(it).estep = estep;
  
      fprintf(['ndfix for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).ndfix = intmp;
      ndfix(it) = intmp;
      
      fprintf(['U for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).U = intmp;
      U(it) = intmp;
  
      fprintf(['V for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).V = intmp;
      V(it) = intmp;
  
      fprintf(['t_0 for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).t_0 = intmp;
      t_0(it) = intmp;

      fprintf(['phi/pi for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');
      partot(it).phi_div_pi = intmp;
      phi_div_pi(it) = intmp;
      
  
      ValidIdx = false;
      while ~ValidIdx
  
        fprintf(['Temperature for job #',sprintf('%d',it),'\n']);
        intmp = input('>>> ');
  
        if isnumeric(intmp) && intmp > 0  
          partot(it).T = intmp;
          T(it) = intmp;
          emin = T(it);
          partot(it).emin = emin;    % Minimum absolute value of frequency grid
          ValidIdx = true;
  
        else
          fprintf('WRN: Invalid Input!\n');
        end
      end

      ocont = getAcont(0, 0, 0, 0, 'emin', emin, 'emax', emax, 'estep', estep);   % logarithmic frequency grid
      initSE = cell(1,2);                   % initial self-energy, even & odd
      initSE{1} =-0.1i*zeros(numel(ocont),3,3);  % initial self-energy, even sector
      initSE{2} = -0.1i*zeros(numel(ocont),2,2);  % initial self-energy, odd sector
      partot(it).initSE = initSE;

      %{}
      JobName = ['U=',sprintf('%.15g',partot(it).U), ...
                '_V=',sprintf('%.15g',partot(it).V), ...
                  '_t_0=',sprintf('%.15g',partot(it).t_0), ...
                    '_phi_div_pi=',sprintf('%.15g',partot(it).phi_div_pi), ...
                      '_T=',sprintf('%.15g',partot(it).T), ...
                        '_ndfix=',sprintf('%.15g',partot(it).ndfix)];  
      %}
  
      %{
      JobName = ['U=',sprintf('%.15g',partot(it).U), ...
                '_V=',sprintf('%.15g',partot(it).V), ...
                  '_t_0=',sprintf('%.15g',partot(it).t_0), ...
                    '_phi_div_pi=',sprintf('%.15g',partot(it).phi_div_pi), ...
                      '_T=',sprintf('%.15g',partot(it).T), ...
                        '_ndfix=',sprintf('%.15g',partot(it).ndfix), ...
                          '_Nfit=27'];     
      %}
  
      partot(it).JobName = JobName;
  
    end
  
    for it = (1:num_jobs)
      if ~exist(['/data/',getenv('USER'),'/DMFT_QI_2CK_lat/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)],'dir')
        mkdir(['/data/',getenv('USER'),'/DMFT_QI_2CK_lat/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)]);
      end
    end
  
    if all(U == U(end)) && num_jobs ~= 1
      Us = ['[',sprintf('%.15g',U(end)),']x',sprintf('%d',num_jobs)];
    else
      Us = cellfun(@(x) [sprintf('%.15g',partot(x).U),','],num2cell(1:num_jobs),'UniformOutput',false);
      Us = cell2mat(Us);
      Us = ['[',Us(1:end-1),']'];
    end
  
    if all(V == V(end)) && num_jobs ~= 1
      Vs = ['[',sprintf('%.15g',V(end)),']x',sprintf('%d',num_jobs)];
    else
      Vs = cellfun(@(x) [sprintf('%.15g',partot(x).V),','],num2cell(1:num_jobs),'UniformOutput',false);
      Vs = cell2mat(Vs);
      Vs = ['[',Vs(1:end-1),']'];
    end
  
    if all(t_0 == t_0(end)) && num_jobs ~= 1
      t_0s = ['[',sprintf('%.15g',t_0(end)),']x',sprintf('%d',num_jobs)];
    else
      t_0s = cellfun(@(x) [sprintf('%.15g',partot(x).t_0),','],num2cell(1:num_jobs),'UniformOutput',false);
      t_0s = cell2mat(t_0s);
      t_0s = ['[',t_0s(1:end-1),']'];
    end
  
    if all(phi_div_pi == phi_div_pi(end)) && num_jobs ~= 1
      phi_div_pis = ['[',sprintf('%.15g',phi_div_pi(end)),']x',sprintf('%d',num_jobs)];
    else
      phi_div_pis = cellfun(@(x) [sprintf('%.15g',partot(x).phi_div_pi),','],num2cell(1:num_jobs),'UniformOutput',false);
      phi_div_pis = cell2mat(phi_div_pis);
      phi_div_pis = ['[',phi_div_pis(1:end-1),']'];
    end
  
    if all(T == T(end)) && num_jobs ~= 1
      Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
    else
      Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
      Ts = cell2mat(Ts);
      Ts = ['[',Ts(1:end-1),']'];
    end

    if all(ndfix == ndfix(end)) && num_jobs ~= 1
      ndfixs = ['[',sprintf('%.15g',ndfix(end)),']x',sprintf('%d',num_jobs)];
    else
      ndfixs = cellfun(@(x) [sprintf('%.15g',partot(x).ndfix),','],num2cell(1:num_jobs),'UniformOutput',false);
      ndfixs = cell2mat(ndfixs);
      ndfixs = ['[',ndfixs(1:end-1),']'];
    end
  
    parfn = ['DMFT_QI_par_U=',Us,'_V=',Vs,'_t_0=',t_0s,'_pdp=',phi_div_pis,'_T=',Ts,'_nu=',ndfixs];
    %parfn = ['DMFT_QI_par_U=',Us,'_V=',Vs,'_t_0=',t_0s,'_pdp=',phi_div_pis,'_T=',Ts,'_nu=',ndfixs,'_Nfit=[27]x2'];
  
    dispstruct(partot);
    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
    disp(['Saved to : ', parfn]);
  end