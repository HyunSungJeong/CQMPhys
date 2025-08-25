clear;

%% (J0, K0, I0) == (-0.08, 0.3, 0)

FitInfo.Imp{1}.Range = [-15,-10];
FitInfo.Imp{1}.LineShift = 4;
FitInfo.Imp{1}.TextShift = [0, 1.8];
FitInfo.Imp{2}.Range = [-15,-10];
FitInfo.Imp{2}.LineShift = 4;
FitInfo.Imp{2}.TextShift = [0, 1.9];
FitInfo.Imp{3}.Range = [-15,-10];
FitInfo.Imp{3}.LineShift = 1/4;
FitInfo.Imp{3}.TextShift = [0, -2.8];

FitInfo.Bath{1}.Range = [-17,-12];
FitInfo.Bath{1}.LineShift = 4;
FitInfo.Bath{1}.TextShift = [0, 2.1];
FitInfo.Bath{2}.Range = [-17,-12];
FitInfo.Bath{2}.LineShift = 4;
FitInfo.Bath{2}.TextShift = [0, 2.3];
FitInfo.Bath{3}.Range = [-10,-6];
FitInfo.Bath{3}.LineShift = 1/4;
FitInfo.Bath{3}.TextShift = [0, -2.8];

YLims = [1e-6, 1e16;
         1e-17, 5e14];

SubFig_idx = {'$\mathrm{(a)}$', '$\mathrm{(b)}$', '$\mathrm{(c)}$'};

alphaz = 1;

NRG_RGflow_fig(-0.08, 0.3, 0, FitInfo, 'YLims', YLims, 'SubFig_idx', SubFig_idx, 'alphaz', alphaz);

%% (J0, K0, I0) == (0.08, 0.3, 0)

FitInfo.Imp{1}.Range = [-18.5,-15;
                         -8.5,-5.5];
FitInfo.Imp{1}.LineShift = [4, 4];
FitInfo.Imp{1}.TextShift = [-0.2, 1.8;
                             0, 1.8];
FitInfo.Imp{2}.Range = [-16,-12];
FitInfo.Imp{2}.LineShift = 4;
FitInfo.Imp{2}.TextShift = [-0.2, 1.8];
FitInfo.Imp{3}.Range = [-17,-14.5];
FitInfo.Imp{3}.LineShift = 1/4;
FitInfo.Imp{3}.TextShift = [0.2, -1.8];

FitInfo.Bath{1}.Range = [-18.5,-15;
                         -8.5,-5.5];
FitInfo.Bath{1}.LineShift = [4, 4];
FitInfo.Bath{1}.TextShift = [-0.2, 2.2;
                             0, 2.2];
FitInfo.Bath{2}.Range = [-16,-12];
FitInfo.Bath{2}.LineShift = 4;
FitInfo.Bath{2}.TextShift = [-0.2, 2.2];
FitInfo.Bath{3}.Range = [-11,-7];
FitInfo.Bath{3}.LineShift = 1/4;
FitInfo.Bath{3}.TextShift = [0, -2.8];

YLims = [1e-6, 1e12;
         1e-16, 1e12];

SubFig_idx = {'$\mathrm{(d)}$', '$\mathrm{(e)}$', '$\mathrm{(f)}$'};

alphaz = 1;

NRG_RGflow_fig(0.08, 0.3, 0, FitInfo, 'YLims', YLims, 'SubFig_idx', SubFig_idx, 'alphaz', alphaz);

%% (J0, K0, I0) == (0.3, 0.08, 0)

FitInfo.Imp{1}.Range = [-16,-12];
FitInfo.Imp{1}.LineShift = 4;
FitInfo.Imp{1}.TextShift = [-0.2, 1.8];
FitInfo.Imp{2}.Range = [-18.5,-15;
                         -8.5,-5.5];
FitInfo.Imp{2}.LineShift = [4, 4];
%FitInfo.Imp{2}.TextShift = [-0.2, 1.8;
%                             0, 1.8];
FitInfo.Imp{2}.TextShift = [-0.2, 1.8;
                             -0.4, 2.2];
FitInfo.Imp{3}.Range = [-17,-14.5];
FitInfo.Imp{3}.LineShift = 1/4;
FitInfo.Imp{3}.TextShift = [0.2, -1.8];

FitInfo.Bath{1}.Range = [-16,-12];
FitInfo.Bath{1}.LineShift = 4;
FitInfo.Bath{1}.TextShift = [-0.2, 2.2];
FitInfo.Bath{2}.Range = [-18.5,-15;
                         -8.5,-5.5];
FitInfo.Bath{2}.LineShift = [4, 4];
%FitInfo.Bath{2}.TextShift = [-0.2, 2.2;
%                             0, 2.2];
FitInfo.Bath{2}.TextShift = [-0.2, 2.2;
                             -0.4, 2.6];
FitInfo.Bath{3}.Range = [-11,-7];
FitInfo.Bath{3}.LineShift = 1/4;
FitInfo.Bath{3}.TextShift = [0, -2.8];

YLims = [1e-6, 1e16;
         1e-16, 5e14];

SubFig_idx = {'$\mathrm{(d)}$', '$\mathrm{(e)}$', '$\mathrm{(f)}$'};

alphaz = 1;

NRG_RGflow_fig(0.3, 0.08, 0, FitInfo, 'YLims', YLims, 'SubFig_idx', SubFig_idx, 'alphaz', alphaz);

%% (J0, K0, I0) == (0.3, -0.08, 0)

FitInfo.Imp{1}.Range = [-15,-10];
FitInfo.Imp{1}.LineShift = 4;
FitInfo.Imp{1}.TextShift = [0, 1.9];
FitInfo.Imp{2}.Range = [-15,-10];
FitInfo.Imp{2}.LineShift = 4;
FitInfo.Imp{2}.TextShift = [0, 1.8];
FitInfo.Imp{3}.Range = [-15,-10];
FitInfo.Imp{3}.LineShift = 1/4;
FitInfo.Imp{3}.TextShift = [0, -2.8];

FitInfo.Bath{1}.Range = [-17,-12];
FitInfo.Bath{1}.LineShift = 4;
FitInfo.Bath{1}.TextShift = [0, 2.3];
FitInfo.Bath{2}.Range = [-17,-12];
FitInfo.Bath{2}.LineShift = 4;
FitInfo.Bath{2}.TextShift = [0, 2.1];
FitInfo.Bath{3}.Range = [-10,-6];
FitInfo.Bath{3}.LineShift = 1/4;
FitInfo.Bath{3}.TextShift = [0, -2.8];

YLims = [1e-6, 1e16;
         1e-17, 5e14];

SubFig_idx = {'$\mathrm{(a)}$', '$\mathrm{(b)}$', '$\mathrm{(c)}$'};

alphaz = 1;

NRG_RGflow_fig(0.3, -0.08, 0, FitInfo, 'YLims', YLims, 'SubFig_idx', SubFig_idx, 'alphaz', alphaz);


%% (J0, K0, I0) == (0.08, 0.3, 1e-6)

FitInfo.Imp{1}.Range = [-18,-14;
                         -8.5, -5];
FitInfo.Imp{1}.LineShift = [4, 4];
FitInfo.Imp{1}.TextShift = [-0.2, 2.2;
                             -0.2, 2.2];
FitInfo.Imp{2}.Range = [-18,-14;
                         -12.4, -10.8;
                         -9.5,-6];
FitInfo.Imp{2}.LineShift = [1/4, 1/4, 4];
FitInfo.Imp{2}.TextShift = [-0.2, -2.7;
                             0, -2.3;
                             -1.6, 2];
FitInfo.Imp{3}.Range = [-18,-14;
                        -8.5, -5];
FitInfo.Imp{3}.LineShift = [4, 1/4];
FitInfo.Imp{3}.TextShift = [-0.2, 2.2;
                            -0.2, -3];

FitInfo.Bath{1}.Range = [-18,-14;
                         -8.5, -5];
FitInfo.Bath{1}.LineShift = [4, 4];
FitInfo.Bath{1}.TextShift = [-0.2, 2;
                             -0.2, 2];
FitInfo.Bath{2}.Range = [-18,-14;
                         -12.4, -10.8;
                         -9.5,-6];
FitInfo.Bath{2}.LineShift = [4, 4, 1/4];
FitInfo.Bath{2}.TextShift = [-0.2, 2;
                             -0.8, 1.8;
                             -0.2, -2.7];
FitInfo.Bath{3}.Range = [-9,-5];
FitInfo.Bath{3}.LineShift = 1/4;
FitInfo.Bath{3}.TextShift = [0, -2.1];

YLims = [1e-16, 1e10;
         1e-16, 1e8];

SubFig_idx = {'$\mathrm{(g)}$', '$\mathrm{(h)}$', '$\mathrm{(i)}$'};

alphaz = [1, 0.8];

NRG_RGflow_fig(0.08, 0.3, 1e-6, FitInfo, 'YLims', YLims, 'SubFig_idx', SubFig_idx, 'alphaz', alphaz);


%% (J0, K0, I0) == (0.3, 0.08, 1e-6)

FitInfo.Imp{1}.Range = [-18,-14;
                         -12.4, -10.8;
                         -9.5,-6];
FitInfo.Imp{1}.LineShift = [1/4, 1/4, 4];
FitInfo.Imp{1}.TextShift = [-0.2, -2.7;
                             0, -2.3;
                             -1.6, 2];
FitInfo.Imp{2}.Range = [-18,-14;
                         -8.5, -5];
FitInfo.Imp{2}.LineShift = [4, 4];
FitInfo.Imp{2}.TextShift = [-0.2, 2.2;
                             -0.2, 2.2];
FitInfo.Imp{3}.Range = [-18,-14;
                        -8.5, -5];
FitInfo.Imp{3}.LineShift = [4, 1/4];
FitInfo.Imp{3}.TextShift = [-0.2, 2.2;
                            -0.2, -3];

FitInfo.Bath{1}.Range = [-18,-14;
                         -12.4, -10.8;
                         -9.5,-6];
FitInfo.Bath{1}.LineShift = [4, 4, 1/4];
FitInfo.Bath{1}.TextShift = [-0.2, 2;
                             -0.8, 1.8;
                             -0.2, -2.7];
FitInfo.Bath{2}.Range = [-18,-14;
                         -8.5, -5];
FitInfo.Bath{2}.LineShift = [4, 4];
FitInfo.Bath{2}.TextShift = [-0.2, 2;
                             -0.2, 2];
FitInfo.Bath{3}.Range = [-9,-5];
FitInfo.Bath{3}.LineShift = 1/4;
FitInfo.Bath{3}.TextShift = [0, -2.1];

YLims = [1e-16, 1e10;
         1e-16, 1e8];

SubFig_idx = {'$\mathrm{(g)}$', '$\mathrm{(h)}$', '$\mathrm{(i)}$'};

alphaz = [1, 0.8];

NRG_RGflow_fig(0.3, 0.08, 1e-6, FitInfo, 'YLims', YLims, 'SubFig_idx', SubFig_idx, 'alphaz', alphaz);


%% (J0, K0, I0) == (0.08, 0.3, 1e-6)

FitInfo.Imp{1}.Range = [-18,-14;
                         -8.5, -5];
FitInfo.Imp{1}.LineShift = [4, 4];
FitInfo.Imp{1}.TextShift = [-0.2, 2.2;
                             -0.2, 2.2];
FitInfo.Imp{2}.Range = [-18,-14;
                         -12.4, -10.8;
                         -9.5,-6];
FitInfo.Imp{2}.LineShift = [1/4, 1/4, 4];
FitInfo.Imp{2}.TextShift = [-0.2, -2.7;
                             0, -2.3;
                             -1.6, 2];
FitInfo.Imp{3}.Range = [-18,-14;
                        -8.5, -5];
FitInfo.Imp{3}.LineShift = [4, 1/4];
FitInfo.Imp{3}.TextShift = [-0.2, 2.2;
                            -0.2, -3];

FitInfo.Bath{1}.Range = [-18,-14;
                         -8.5, -5];
FitInfo.Bath{1}.LineShift = [4, 4];
FitInfo.Bath{1}.TextShift = [-0.2, 2;
                             -0.2, 2];
FitInfo.Bath{2}.Range = [-18,-14;
                         -12.4, -10.8;
                         -9.5,-6];
FitInfo.Bath{2}.LineShift = [4, 4, 1/4];
FitInfo.Bath{2}.TextShift = [-0.2, 2;
                             -0.8, 1.8;
                             -0.2, -2.7];
FitInfo.Bath{3}.Range = [-9,-5];
FitInfo.Bath{3}.LineShift = 1/4;
FitInfo.Bath{3}.TextShift = [0, -2.1];

YLims = [1e-18, 1e12;
         1e-18, 1e12];

SubFig_idx = {'$\mathrm{(g)}$', '$\mathrm{(h)}$', '$\mathrm{(i)}$'};

alphaz = [1, 0.8];

NRG_RGflow_fig(0.08, 0.3, 1e-6, FitInfo, 'YLims', YLims, 'SubFig_idx', SubFig_idx, 'alphaz', alphaz);