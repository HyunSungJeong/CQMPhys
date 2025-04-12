function TsoK_Aniso_par (varargin)

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
  h_vmem = 64;           % Memory (in GB) to be occupied in clusters
  PE = 7;                 % # of cores to be occupied in clusters
  syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
  Nkeep = 3000;
  nz = ones(1,num_jobs);
  J0 = zeros(1,num_jobs);       % spin exchange coupling
  K_perp = zeros(1,num_jobs);   % orbital pseudospin exchange coupling: perpendicular component
  K_z = zeros(1,num_jobs);      % orbital pseudospin exchange coupling: z-component
  I0 = zeros(1,num_jobs);       % spin-orbital exchange coupling
  T = zeros(1,num_jobs);      % temperature

  for it = (1:num_jobs)

    partot(it).nz = nz(it);
    partot(it).PE = PE;
    partot(it).Nkeep = Nkeep;

    fprintf(['J0 for job #',sprintf('%d',it),'\n']);
    intmp = input('>>> ');
    partot(it).J0 = intmp;
    J0(it) = intmp;
    
    fprintf(['K_perp for job #',sprintf('%d',it),'\n']);
    intmp = input('>>> ');
    partot(it).K_perp = intmp;
    K_perp(it) = intmp;

    fprintf(['K_z for job #',sprintf('%d',it),'\n']);
    intmp = input('>>> ');
    partot(it).K_z = intmp;
    K_z(it) = intmp;

    fprintf(['I0 for job #',sprintf('%d',it),'\n']);
    intmp = input('>>> ');
    partot(it).I0 = intmp;
    I0(it) = intmp;

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

    JobName = ['J0=',sprintf('%.15g',partot(it).J0), ...
              '_K_perp=',sprintf('%.15g',partot(it).K_perp), ...
                '_K_z=',sprintf('%.15g',partot(it).K_z), ...
                  '_I0=',sprintf('%.15g',partot(it).I0), ...
                    '_T=',sprintf('%.15g',partot(it).T)];     

    partot(it).JobName = JobName;

  end

  for it = (1:num_jobs)
    if ~exist(['/data/',getenv('USER'),'/TsoK_Aniso/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)],'dir')
      mkdir(['/data/',getenv('USER'),'/TsoK_Aniso/',partot(it).JobName,'_Nkeep=',sprintf('%.15g',Nkeep)]);
    end
  end

  if all(J0 == J0(end)) && num_jobs ~= 1
    J0s = ['[',sprintf('%.15g',J0(end)),']x',sprintf('%d',num_jobs)];
  else
    J0s = cellfun(@(x) [sprintf('%.15g',partot(x).J0),','],num2cell(1:num_jobs),'UniformOutput',false);
    J0s = cell2mat(J0s);
    J0s = ['[',J0s(1:end-1),']'];
  end

  if all(K_perp == K_perp(end)) && num_jobs ~= 1
    K_perps = ['[',sprintf('%.15g',K_perp(end)),']x',sprintf('%d',num_jobs)];
  else
    K_perps = cellfun(@(x) [sprintf('%.15g',partot(x).K_perp),','],num2cell(1:num_jobs),'UniformOutput',false);
    K_perps = cell2mat(K_perps);
    K_perps = ['[',K_perps(1:end-1),']'];
  end

  if all(K_z == K_z(end)) && num_jobs ~= 1
    K_zs = ['[',sprintf('%.15g',K_z(end)),']x',sprintf('%d',num_jobs)];
  else
    K_zs = cellfun(@(x) [sprintf('%.15g',partot(x).K_z),','],num2cell(1:num_jobs),'UniformOutput',false);
    K_zs = cell2mat(K_zs);
    K_zs = ['[',K_zs(1:end-1),']'];
  end

  if all(I0 == I0(end)) && num_jobs ~= 1
    I0s = ['[',sprintf('%.15g',I0(end)),']x',sprintf('%d',num_jobs)];
  else
    I0s = cellfun(@(x) [sprintf('%.15g',partot(x).I0),','],num2cell(1:num_jobs),'UniformOutput',false);
    I0s = cell2mat(I0s);
    I0s = ['[',I0s(1:end-1),']'];
  end

  if all(T == T(end)) && num_jobs ~= 1
    Ts = ['[',sprintf('%.15g',T(end)),']x',sprintf('%d',num_jobs)];
  else
    Ts = cellfun(@(x) [sprintf('%.15g',partot(x).T),','],num2cell(1:num_jobs),'UniformOutput',false);
    Ts = cell2mat(Ts);
    Ts = ['[',Ts(1:end-1),']'];
  end

  parfn = ['TsoK_Aniso_par_J0=',J0s,'_K_perp=',K_perps,'_K_z=',K_zs,'_I0=',I0s,'_T=',Ts];

  dispstruct(partot);
  parfn = [go('mu/Para/'), parfn, '.mat'];
  save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
  disp(['Saved to : ', parfn]);
end