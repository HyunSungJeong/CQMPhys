function bathNRG_par (varargin)

  while numel(varargin) > 0

    if ischar(varargin{1})
      parfn = varargin{1};
      varargin(1) = [];
    else
      disp(varargin{1}); 
      error('ERR: Unknown input/option.');
    end

  end

  syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
  h_vmem = 250;           % Memory (in GB) to be occupied in clusters
  PE = 20;                 % # of cores to be occupied in clusters
  syms = cell(1, 0);      % non-Abelian symmetry types to be exploited
  T = 1e-30;              % Temperature

  parfn = ['bathNRG_par_T=',sprintf('%.15g',T)];

  partot.PE = PE;
  partot.T = T;

  dispstruct(partot);
  parfn = [go('mu/Para/'), parfn, '.mat'];
  %parfn = ['/home/',getenv('USER'),'/Para/',parfn,'.mat'];
  save(parfn, 'partot', 'PE', 'h_vmem', 'syms');
  disp(['Saved to : ', parfn]);
end