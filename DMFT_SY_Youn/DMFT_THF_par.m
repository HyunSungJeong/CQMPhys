function DMFT_THF_par (varargin)

    parfn = 'DMFT_AHM_par'; % default file name

    while numel(varargin) > 0

        if ischar(varargin{1})
        parfn = varargin{1};
        varargin(1) = [];
        else
        disp(varargin{1});
        error('ERR: Unknown input/option.');
        end

    end

    h_vmem = 200; % Memory (in GB) to be occupied in clusters
    PE = 10; % # of cores to be occupied in clusters
    syms = cell(1, 0); % non-Abelian symmetry types to be exploited

    U = [1:0.5:6];
    % U = 1;
    delta = 0.2;
    r = 0.3;
    % muc = U / 2;
    ndfix = 2; %total filling = 0 from CNP, half filling
    T = [1, 1e-2, 1e-4];
    Nkeep = 1200;

    initS = [];
    % %% initfile scheme
    % initfile = [
    %             "/home/seongyeon/AHM/AHM_T=0.0001_U=1_r=0.3_Lambda=3_nz=2_Nk=1200_j13929t1.mat";
    %             "/home/seongyeon/AHM/AHM_T=0.0001_U=2_r=0.3_Lambda=3_nz=2_Nk=1200_j13929t2.mat";
    %             "/home/seongyeon/AHM/AHM_T=0.0001_U=2.4_r=0.3_Lambda=3_nz=2_Nk=1200_j13929t3.mat";
    %             "/home/seongyeon/AHM/AHM_T=0.0001_U=2.8_r=0.3_Lambda=3_nz=2_Nk=1200_j13929t4.mat";
    %             "/home/seongyeon/AHM/AHM_T=0.0001_U=3.2_r=0.3_Lambda=3_nz=2_Nk=1200_j13929t5.mat";
    %             "/home/seongyeon/AHM/AHM_T=0.0001_U=4_r=0.3_Lambda=3_nz=2_Nk=1200_j13929t6.mat";
    %             "/home/seongyeon/AHM/AHM_T=0.0001_U=6_r=0.3_Lambda=3_nz=2_Nk=1200_j13929t7.mat";
    %             ];

    % initS = cell(numel(initfile), 1);

    % for l = (1:numel(initfile))
    %   initS{l} = struct;
    %   D = load(initfile(l), 'SEs', 'it', 'mures', 'ocont');
    %   initS{l}.SE{1} = D.SEs{1}(:, :, D.it);
    %   initS{l}.SE{2} = D.SEs{2}(:, :, D.it);
    %   initS{l}.mu = D.mures(D.it);
    %   initS{l}.ocont = D.ocont;
    % end

    Lambda = 3;
    nz = 2;
    sigmarat = 1.5;

    partot = struct;
    a = 0;

    for i = (1:numel(U))

        for j = (1:numel(r))

        for k = (1:numel(T))

            a = a + 1;
            partot(a).PE = PE;
            % partot(a).muc = muc(i);
            partot(a).ndfix = ndfix;
            partot(a).T = T(k);
            partot(a).Lambda = Lambda;
            partot(a).nz = nz;
            partot(a).sigmarat = sigmarat;
            partot(a).Nkeep = Nkeep;
            partot(a).U = U(i);
            partot(a).r = r(j);
            partot(a).delta = delta;
            partot(a).initS = [];
            % partot(a).initS = initS{i};
        end

        end

    end

    partot = partot(:);

    dispstruct(partot);

    disp(partot(1))

    parfn = [go('mu/Para/'), parfn, '.mat'];
    save(parfn, 'partot', 'h_vmem', 'PE', 'syms', '-v7.3');
    disp(['Saved to : ', parfn]);
end