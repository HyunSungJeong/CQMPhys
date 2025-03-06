function LineSearch_par (varargin)

  syms = cell(1, 0);            % non-Abelian symmetry types to be exploited
  h_vmem = 250;                 % Memory (in GB) to be occupied in clusters
  PE = 10;                      % # of cores to be occupied in clusters
  syms = cell(1, 0);            % non-Abelian symmetry types to be exploited
  T = 1e-24;

  fprintf('Number of Jobs?\n');
  num_jobs = input('>>> ');

  for it = (1:num_jobs)
    partot(it).PE = PE;
    partot(it).T = T;
    partot(it).Bi_Par = [];
    patot(it).J0 = [];
    partot(it).K0 = [];
    partot(it).I0 = [];
    partot(it).hi = [];
    partot(it).lo = [];
  end

  Bi_Par = [];
  hi = [];
  lo = [];

  while isempty(Bi_Par)
    fprintf('Index for parameter to be binary-searched:\n');
    fprintf('1: J0\n');
    fprintf('2: K0\n');
    fprintf('3: I0\n');
    intmp = input('>>> ');
    
    if intmp == 1
      Bi_Par = 'J0';
    elseif intmp == 2
      Bi_Par = 'K0';
    elseif intmp == 3
      Bi_Par = 'I0';
    else
      fprintf('WRN: invalid input\n');
    end


  for it = (1:num_jobs)
    while isempty(partot(it).Bi_Par)

      if Bi_Par == 'J0'

        partot(it).Bi_Par = 'J0';
        fprintf(['Fixed value of K0 for job #',sprintf('%d',it),' :\n']);
        intmp = input('>>> ');
        partot(it).K0 = intmp;
        fprintf(['Fixed value of I0 for job #',sprintf('%d',it),' :\n']);
        intmp = input('>>> ');
        partot(it).I0 = intmp;

      elseif Bi_Par == 'K0'

        partot(it).Bi_Par = 'K0';
        fprintf(['Fixed value of J0 for job #',sprintf('%d',it),' :\n']);
        intmp = input('>>> ');
        partot(it).J0 = intmp;
        fprintf(['Fixed value of I0 for job #',sprintf('%d',it),' :\n']);
        intmp = input('>>> ');
        partot(it).I0 = intmp;

      elseif Bi_Par == 'I0'

        partot(it).Bi_Par = 'I0';
        fprintf(['Fixed value of J0 for job #',sprintf('%d',it),' :\n']);
        intmp = input('>>> ');
        partot(it).J0 = intmp;
        fprintf(['Fixed value of K0 for job #',sprintf('%d',it),' :\n']);
        intmp = input('>>> ');
        partot(it).K0 = intmp;

      end
    end
  end

  while isempty(hi)
    fprintf(['Upper bound for ',Bi_Par,':\n']);
    intmp1 = input('>>> ');
    fprintf(['Lower bound for ',Bi_Par,':\n']);
    intmp2 = input('>>> ');

    if intmp1 > intmp2
      hi = intmp1;
      lo = intmp2;
    else
      fprintf('WRN: Upper bound must be greater than lower bound\n');
    end
  end

  for it = (1:num_jobs)
    partot(it).hi = hi;
    partot(it).lo = lo;
  end

  if num_jobs == 1

    if Bi_Par == 'J0'

      partot(1).JobName = ['J0_search_',sprintf('%.15g',partot(1).lo),'_to_',sprintf('%.15g',partot(1).hi),'_K0=',sprintf('%.15g',partot(1).K0),'_I0=',sprintf('%.15g',partot(1).I0)];
      partot(1).JobName = [partot.JobName,'_MaxSep=1e-2_NumComp=4'];    

    elseif Bi_Par == 'K0'

      partot(1).JobName = ['K0_search_',sprintf('%.15g',partot(1).lo),'_to_',sprintf('%.15g',partot(1).hi),'_J0=',sprintf('%.15g',partot(1).J0),'_I0=',sprintf('%.15g',partot(1).I0)];

    elseif Bi_Par == 'I0'

      partot(1).JobName = ['I0_search_',sprintf('%.15g',partot(1).lo),'_to_',sprintf('%.15g',partot(1).hi),'_J0=',sprintf('%.15g',partot(1).J0),'_K0=',sprintf('%.15g',partot(1).K0)];

    end

    parfn = ['LineSearch_par_',partot(1).JobName];
  else
    for it = (1:num_jobs)

      if Bi_Par == 'J0'

        partot(it).JobName = ['J0_search_',sprintf('%.15g',lo),'_to_',sprintf('%.15g',hi),'_K0=',sprintf('%.15g',partot(it).K0),'_I0=',sprintf('%.15g',partot(it).I0)];
        partot(it).JobName = [partot.JobName,'_MaxSep=1e-2_NumComp=4'];    

      elseif Bi_Par == 'K0'

        partot(it).JobName = ['K0_search_',sprintf('%.15g',lo),'_to_',sprintf('%.15g',hi),'_J0=',sprintf('%.15g',partot(it).J0),'_I0=',sprintf('%.15g',partot(it).I0)];

      elseif Bi_Par == 'I0'

        partot(it).JobName = ['I0_search_',sprintf('%.15g',lo),'_to_',sprintf('%.15g',hi),'_J0=',sprintf('%.15g',partot(it).J0),'_K0=',sprintf('%.15g',partot(it).K0)];

      end

    end

    if Bi_Par == 'J0'  

      K0s = cellfun(@(x) [sprintf('%.15g',partot(x).K0),','],num2cell(1:num_jobs),'UniformOutput',false);
      K0s = cell2mat(K0s);
      K0s = K0s(1:end-1)
      I0s = cellfun(@(x) [sprintf('%.15g',partot(x).I0),','],num2cell(1:num_jobs),'UniformOutput',false);
      I0s = cell2mat(I0s);
      I0s = I0s(1:end-1);
      parfn = ['LineSearch_par_J0_search_',sprintf('%.15g',lo),'_to_',sprintf('%.15g',hi), ...
                  '_K0=[',K0s,']_I0=[',I0s,']'];

    elseif Bi_Par == 'K0'

      J0s = cellfun(@(x) [sprintf('%.15g',partot(x).J0),','],num2cell(1:num_jobs),'UniformOutput',false);
      J0s = cell2mat(J0s);
      J0s = J0s(1:end-1)
      I0s = cellfun(@(x) [sprintf('%.15g',partot(x).I0),','],num2cell(1:num_jobs),'UniformOutput',false);
      I0s = cell2mat(I0s);
      I0s = I0s(1:end-1);
      parfn = ['LineSearch_par_K0_search_',sprintf('%.15g',lo),'_to_',sprintf('%.15g',hi), ...
                  '_J0=[',J0s,']_I0=[',I0s,']'];

    elseif Bi_Par == 'I0'

      J0s = cellfun(@(x) [sprintf('%.15g',partot(x).J0),','],num2cell(1:num_jobs),'UniformOutput',false);
      J0s = cell2mat(J0s);
      J0s = J0s(1:end-1)
      K0s = cellfun(@(x) [sprintf('%.15g',partot(x).K0),','],num2cell(1:num_jobs),'UniformOutput',false);
      K0s = cell2mat(K0s);
      K0s = K0s(1:end-1);
      parfn = ['LineSearch_par_I0_search_',sprintf('%.15g',lo),'_to_',sprintf('%.15g',hi), ...
                  '_J0=[',J0s,']_K0=[',K0s,']'];

    end
  end

  for it = (1:num_jobs)
    if ~exist(['/data/',getenv('USER'),'/LineSearch/',Bi_Par,'_search/',partot(it).JobName],'dir')
          mkdir(['/data/',getenv('USER'),'/LineSearch/',Bi_Par,'_search/',partot(it).JobName]);
    end
  end

  dispstruct(partot);

    while numel(varargin) > 0

        if ischar(varargin{1})
        parfn = varargin{1};
        varargin(1) = [];
        else
        disp(varargin{1}); 
        error('ERR: Unknown input/option.');
        end
    end

  parfn = [go('mu/Para/'), parfn, '.mat'];
  %parfn = ['/home/',getenv('USER'),'/Para/',parfn,'.mat'];
  save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
  disp(['Saved to : ', parfn]);
end