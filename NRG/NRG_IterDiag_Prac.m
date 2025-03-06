function Inrg = NRG_IterDiag(H0,A0,Lambda,ff,F,gg,NF,Z,Nkeep)
    
    Nfac = 0.1;

    % % error checking
    if numel(H0.Q) ~= 2
        error('ERR: ''H0'' should be of rank 2.');
    elseif numel(A0.Q) ~= 3
        error('ERR: ''A0'' should be of rank 3.');
    end

    Inrg = struct;
    Inrg.Lambda = Lambda;
    N = numel(ff) + 1;

    Inrg.EScale = [1, (Lambda.^(((N-2):-1:0)/2))*ff(end)];

    % NRG results
    Inrg.EK = cell(1,N);
    Inrg.AK = cell(1,N);
    Inrg.ED = cell(1,N);
    Inrg.AD = cell(1,N);
    Inrg.E0 = zeros(1,N);

    tobj = tic2;

    disptime('NRG: start');

    for itN = (1:N)
        if itN == 1
            A0.info.itags = {'imp','H0*','L0'};
            H0.info.itags = {'H0','H0*'};
            Z.info.itags = {'L0','L0*'};
            F.info.itags = {'L0','L0*','op*'};
            NF.info.itags = {'L0','L0*'};

            Inrg.AK{itN} = A0;
            [E,I] = eigQS(H0);
            Ntr = size(E,1);
            Inrg.EK{itN} = QSpace(I.EK);
            Inrg.AD{itN} = QSpace;   % How should I handle these??
            Inrg.ED{itN} = QSpace;

            Hprev = diag(Inrg.EK{itN});
        else
            
            Fprev = contract(F,'!3*',{Inrg.AK{itN-1},'!2*',Inrg.AK{itN-1},'!3'},'!1',[2,3,1]);      % Mind the order of contraction!(zipper)

            Z.info.itags = {append('L',int2str(itN-1)), append('L',int2str(itN-1),'*')};
            F.info.itags = {append('L',int2str(itN-1)), append('L',int2str(itN-1),'*'),'op*'};
            NF.info.itags = {append('L',int2str(itN-1)), append('L',int2str(itN-1),'*')};

            Anow = getIdentity(Inrg.AK{itN-1},2,Z,2,append('H',int2str(itN-1),'*'),[1,3,2]);
            Hnow = contract({Anow,'!2*',{Hprev,'!1',Anow}});
            Hnow = Hnow*Inrg.EScale(itN-1)/Inrg.EScale(itN);

            Hhop = contract(Anow,'!2*',{{Z,'!1',F},{Fprev,Anow}});      % Mind the order of contraction!(zipper)
            Hhop = Hhop*ff(itN-1)/Inrg.EScale(itN);
            Hhop = Hhop + permute(Hhop,[2,1],'conj');

            Hon = contract(NF,{Anow,'!2*',Anow,'!3'},'!1');
            Hon = Hon*gg(itN-1)/Inrg.EScale(itN-1);

            Hnow = Hnow + Hhop + Hon;

            [E,Ieig] = eigQS((Hnow + permute(Hnow,[2,1],'conj'))/2,'Nkeep',Nkeep);%,'deps',10);
            
            Inrg.EK{itN} = QSpace(Ieig.EK);
            Inrg.AK{itN} = QSpace(Ieig.AK);
            Inrg.ED{itN} = QSpace(Ieig.ED);
            Inrg.AD{itN} = QSpace(Ieig.AD);

            
            for it = (1:numel(Inrg.EK{itN}.data))
                Inrg.EK{itN}.data{it} = Inrg.EK{itN}.data{it} - E(1,1);
            end

            for it = (1:numel(Inrg.ED{itN}.data))
                Inrg.ED{itN}.data{it} = Inrg.ED{itN}.data{it} - E(1,1);
            end

            Inrg.AD{itN} = contract(Anow,Inrg.AD{itN},[1,3,2]);
            Inrg.AK{itN} = contract(Anow,Inrg.AK{itN},[1,3,2]);

            Hprev = diag(Inrg.EK{itN});
        end

        % information on truncation
        if isempty(Inrg.EK{itN})
            Etr1 = 0;
            Ntr1 = 0;
        else
            Ntr1 = 0;
            for it = (1:numel(Inrg.EK{itN}.data))
                Ntr1 = Ntr1 + numel(Inrg.EK{itN}.data{it});
            end
            Etr1 = E(Ntr1,1) - E(1,1);
        end
        
        if isempty(Inrg.ED{itN})
            Etr2 = Etr1;
        else
            Etr2 = E(end,1) - E(1,1);
        end
        Ntr2 = size(E,1);
        disptime(['#',sprintf('%02i/%02i',[itN-1,N-1]),' : ', ...
            'NK=',sprintf('%i/%i',[Ntr1,Ntr2]),', ', ...
            'EK=',sprintf('%.4g/%.4g',[Etr1,Etr2])]);
    end

    chkmem;
    toc2(tobj,'-v');
end