function HundKM (parfn, varargin)

    try  % % % in case of bug
    
    
        
    partot = job_func_preamble(parfn, varargin{:});
    
    % START : MAIN BODY OF FUNCTION
    
    for it = (1:numel(partot))
        % % load input parameters
        [ J , K , I , mu , D , T , Lambda , nz,  sigmarat , Nkeep , Etrunc , isBroyden  ,  emin, sym,   etrunc,    nd] = loadvar(partot(it), ...
        {'J','K','I','mu','D', 'T','Lambda','nz','sigmarat','Nkeep','Etrunc','isBroyden','emin','sym','etrunc',   'nd'}, ...
        { [], [], [], [] , [] , [],   []  ,  [] ,    []    ,   []  ,   []   ,    []     ,   [] , []  ,   []  ,     [] });
    
    
        % % other parameters
    
        N = max(ceil(-2*log(T/100)/log(Lambda))+8,10); % have to know 
    
    
    % % result file path
        resname = go(['glo/2orbHundKM/2orbHundKM_',par2str('sym',sym,'J',J,'K',K,'I',I,'T',T,'nd',nd,'Nk',Nkeep,'Lambda',Lambda,'nz',nz)]);
        if ~isempty(getenv('SLURM_JOB_ID'))
            resname = [resname,'_j',getenv('SLURM_ARRAY_JOB_ID'),'t',getenv('SLURM_ARRAY_TASK_ID'),'.mat'];
        else
            resname = [resname,'_',char(datetime('now','Format','yy.MM.dd_HH-mm-ss')),'.mat'];
        end
    
        nrgdata = cellfun(@(x) go(['data/NRG_itz=',sprintf('%i',x)]), num2cell(1:nz), 'UniformOutput', false);
    
    
        % % Define QSpace objects
        if isequal(sym,'su2') % % % % U(1) charge * SU(2) spin * SU(2) orbital
            [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge, SU2spin, SU2channel');
            LF = quadOp(FF,FF,[0 0 2]); % orbital opertors;
            SLF = quadOp(FF,FF,[0 2 2]);
            [Fs,Zs,Ss,Ls, SLs Is] = setItag('s00','op',FF,ZF,SF,LF,SLF, IF.E);
            [ZL,SL,LL,SLL, IL] = setItag('L00','op',ZF,SF,LF,SLF, IF.E);
    
    %         project onto half-filled subspace
            IL = getsub(IL,find(IL.Q{1}(:,1)==0));
            Z_L00 = getsub(ZL,find(ZL.Q{1}(:,1)==0));
            A0 = getIdentity(IL,2,Is,2,'K00*',[1 3 2]); % isometry
    
            HSS = contract(SL, '!2*', Ss, [2 1 3 4]); % spin-spin interaction
            HLL = contract(LL, '!2*', Ls, [2 1 3 4]); % orbital-orbital interaction
            HSL = contract(SLL, '!2*', SLs, [2 1 3 4]); %spin-orb - spin-orb interaction 
    
    %         Htot = HSS+HLL+HSL;
            Htot = J*HSS+K*HLL+I*HSL
    
            Op = [(contract(Fs,'!1',Htot,[3 4 1 5 2]) - contract(Htot,'!13',Fs)); SL; LL; SLL]
            cflag = [1; -1;-1;-1];
            zflag = [1;0;0;0];
    
            H0 = contract(A0,'!2*',{Htot,'!13',A0}) + 1e-40*getIdentity(A0,2);
           
    
        end
    
    
        if isequal(sym,'su3') % % % % U(1) charge * SU(2) spin * SU(2) orbital
            [FF, ZF, SF, IF] = getLocalSpace('FermionS', 'Acharge, SU2spin, SU3channel');
            
    %         zf = getsub(ZF, find(ZF.Q{1}(:,3)==0));
    %         zf = getsub(zf, find(zf.Q{1}(:,4)==1));
            
            sf = getsub(SF,find(SF.Q{1}(:,3)==0)); % spin operators
            sf = getsub(sf,find(sf.Q{1}(:,4)==1));
            sf = getsub(sf,find(sf.Q{2}(:,3)==0));
            sf = getsub(sf,find(sf.Q{2}(:,4)==1));
            
            LF = quadOp(FF,FF,[0 0 1 1]); % orbital opertors;
            lf = getsub(LF,find(LF.Q{1}(:,3)==0));
            lf = getsub(lf,find(lf.Q{1}(:,4)==1));
            lf = getsub(lf,find(lf.Q{2}(:,3)==0));
            lf = getsub(lf,find(lf.Q{2}(:,4)==1));
    
            SLF = quadOp(FF,FF,[0 2 1 1]); % spin-orbital opertors;
            slf = getsub(SLF,find(SLF.Q{1}(:,3)==0));
            slf = getsub(slf,find(slf.Q{1}(:,4)==1));
            slf = getsub(slf,find(slf.Q{2}(:,3)==0));
            slf = getsub(slf,find(slf.Q{2}(:,4)==1));
    
    
            ef = getsub(IF.E,find(IF.E.Q{1}(:,3)==0));
            ef = getsub(ef,find(ef.Q{1}(:,4)==1));
            ef = getsub(ef,find(ef.Q{2}(:,3)==0));
            ef = getsub(ef,find(ef.Q{2}(:,4)==1));
    
            
            [Fs,Zs,Ss,Ls, SLs, Is] = setItag('s00','op',FF,zf,sf,lf,slf,ef);
            [ZL,SL,LL,IL] = setItag('L00','op',ZF,sf,LF, IF.E);
            
    
    
    %         project onto half-1 subspace
            A0 = getIdentity(Is,2,IL,2,'K00*',[1 3 2]); % isometry
    
            HSS = contract(SL, '!2*', Ss, [2 1 3 4]);
            HLL = contract(LL, '!2*', Ls, [2 1 3 4]);
            HSL = contract(SLL, '!2*', SLs, [2 1 3 4]);
    
    
    %         Htot = HSS+HLL+HSL;
            Htot = J*HSS+K*HLL+I*HSL
    
            Op = [(contract(Fs,'!1',Htot,[3 4 1 5 2]) - contract(Htot,'!13',Fs)); SL; LL; SLL];
            cflag = [1; -1;-1;-1];
            zflag = [1;0;0;0];
    
            H0 = contract(A0,'!2*',{Htot,'!13',A0}) + 1e-40*getIdentity(A0,2);
            
    
        end
        
    
    
    ff = doZLD([-1 1],[1 1],Lambda,N,nz);
    
    
    Adisc = cell(1,1,nz,numel(Op));
    Acont = cell(1,numel(Op));
    Aw0 = zeros(1,numel(Op),nz); % static susceptibilities
    
    
    for itz = (1:nz)
        nrgdata = NRG_SL([],H0,A0,Lambda,ff{itz}(2:end),FF,ZF,'Nkeep',Nkeep);
        if itz == nz
            [Es,Qs] = plotE(nrgdata);
        end
        
        nrgdata = getRhoFDM(nrgdata,T,'-v');
        [odisc,Adiscz,sigmak,Aw0(1,:,itz)] = getAdisc(nrgdata,Op(:),Op(:),ZF,'Z_L00',Z_L00,'cflag',cflag,'zflag',zflag);
        if numel(Op) > 1
            Adisc(:,:,itz,:) = Adiscz;
        else
            Adisc(:,:,itz,:) = {Adiscz};
        end
    end
    
    
    % z average
    for ito = (1:numel(Op))
        [ocont,Acont{1,ito}] = getAcont(odisc,mean(cell2mat(Adisc(:,:,:,ito)),3),sigmak,T/5,'alphaz',1/nz);
    end
    
    % figure;
    % plot(ocont(ocont>0),Acont{1}(ocont>0)*(pi^2/2/2)); % /2 due to the sum over spin
    % set(gca,'XScale','log');
    % grid on;
    % xlabel('\omega');
    % ylabel('T(\omega)');
    % title('Spin-1/2 + 1-channel: U(1) charge * SU(2) spin');
    % 
    % figure;
    % plot(ocont(ocont>0),Acont{2}(ocont>0));
    % set(gca,'XScale','log','YScale','log');
    % grid on;
    % xlabel('\omega');
    % ylabel('(-1/\pi) Im \chi (\omega)');
        if ~exist(fileparts(resname),'dir')
                mkdir(fileparts(resname));
        end
        save(resname,'-v7.3');
        disptime(['Saved data to: ',resname]);
    end
    catch e
        disp(getReport(e)); % Report error
        if ~(ismcc || isdeployed)
            keyboard;
        end
        
        % Save current workspace
        errf = ['Error_',mfilename];
        if ~isempty(getenv('SLURM_ARRAY_JOB_ID'))
            errf = [errf,'_j',getenv('SLURM_ARRAY_JOB_ID')];
            if ~isempty(getenv('SLURM_ARRAY_TASK_ID'))
                errf = [errf,'t',getenv('SLURM_ARRAY_TASK_ID')];
            end
        end
        
        errf = [go('glo'),filesep,errf,'.mat'];
        save(errf,'-v7.3');
        dispbox('Error! Current workspace is saved in :',errf);
        
        % everything should be before rethrowing error
        rethrow(e);
    end
    
    
    end