function TsoK_par (varargin)

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
  Nkeep = 3000;
  Lambda = 2.5;           % NRG discretization parameter
  getSusc = true;        % whether to compute dynamic susceptibility or not
  getCorr = false;        % whether to compute correlation functions or not
  nz = ones(1,num_jobs);
  J0 = zeros(1,num_jobs);
  K0 = zeros(1,num_jobs);
  I0 = zeros(1,num_jobs);
  T = zeros(1,num_jobs);
  
  for it = (1:num_jobs)

    partot(it).nz = nz(it);
    partot(it).PE = PE;
    partot(it).Nkeep = Nkeep;
    partot(it).getSusc = getSusc;
    partot(it).getCorr = getCorr;
    partot(it).Lambda = Lambda;

    fprintf(['J0 for job #',sprintf('%d',it),'\n']);
    intmp = input('>>> ');
    partot(it).J0 = intmp;
    J0(it) = intmp;
    
    fprintf(['K0 for job #',sprintf('%d',it),'\n']);
    intmp = input('>>> ');
    partot(it).K0 = intmp;
    K0(it) = intmp;

    fprintf(['I0 for job #',sprintf('%d',it),'\n']);
    intmp = input('>>> ');
    partot(it).I0 = intmp;
    I0(it) = intmp;

    ValidIdx = false;
    while ~ValidIdx

      fprintf(['Temperature for job #',sprintf('%d',it),'\n']);
      intmp = input('>>> ');

      if isnumeric(intmp) && intmp > 1e-45 
        partot(it).T = intmp;
        ValidIdx = true;
        T(it) = intmp;
      else
        fprintf('WRN: Invalid Input!\n');
      end
    end

    JobName = ['J0=',sprintf('%.15g',partot(it).J0), ...
              '_K0=',sprintf('%.15g',partot(it).K0), ...
                '_I0=',sprintf('%.15g',partot(it).I0), ...
                  '_T=',sprintf('%.15g',partot(it).T)];     

    partot(it).JobName = JobName;

  end

  for it = (1:num_jobs)
    if ~exist(['/data/',getenv('USER'),'/TsoK/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)],'dir')

      if exist(['/data/',getenv('USER'),'/TsoK/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)],'dir')

        movefile(['/data/',getenv('USER'),'/TsoK/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)], ...
                  ['/data/',getenv('USER'),'/TsoK/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)]);
      else
        mkdir(['/data/',getenv('USER'),'/TsoK/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)]);
      end
    end
  end

  if all(J0 == J0(end))
    J0s = ['[',sprintf('%.15g',J0(end)),']x',sprintf('%d',num_jobs)];
  else
    J0s = cellfun(@(x) [sprintf('%.15g',partot(x).J0),','],num2cell(1:num_jobs),'UniformOutput',false);
    J0s = cell2mat(J0s);
    J0s = ['[',J0s(1:end-1),']'];
  end

  if all(K0 == K0(end))
    K0s = ['[',sprintf('%.15g',K0(end)),']x',sprintf('%d',num_jobs)];
  else
    K0s = cellfun(@(x) [sprintf('%.15g',partot(x).K0),','],num2cell(1:num_jobs),'UniformOutput',false);
    K0s = cell2mat(K0s);
    K0s = ['[',K0s(1:end-1),']'];
  end

  if all(I0 == I0(end))
    I0s = ['[',sprintf('%.15g',I0(end)),']x',sprintf('%d',num_jobs)];
  else
    I0s = cellfun(@(x) [sprintf('%.15g',partot(x).I0),','],num2cell(1:num_jobs),'UniformOutput',false);
    I0s = cell2mat(I0s);
    I0s = ['[',I0s(1:end-1),']'];
  end

  if all(T == T(end))
    Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
  else
    Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
    Ts = cell2mat(Ts);
    Ts = ['[',Ts(1:end-1),']'];
  end

  parfn = ['TsoK_par_J0=',J0s,'_K0=',K0s,'_I0=',I0s,'_T=',Ts,'_Nkeep=',sprintf('%.15g',Nkeep),'_Lambda=',sprintf('%.15g',Lambda)];

  dispstruct(partot);
  parfn = [go('mu/Para/'), parfn, '.mat'];
  save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
  disp(['Saved to : ', parfn]);
end