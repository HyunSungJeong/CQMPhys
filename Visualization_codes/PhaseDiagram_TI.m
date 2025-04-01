function PhaseDiagram_TI(J0,K0,varargin)
    % <Description>
    % Plots T-I phase diagram of the isotropic 2soK model for given (J0,K0)
    %
    % <Input>
    % J0 : [numeric] spin-spin coupling for the isotropic 2-channel spin-orbital Kondo(2soK) model 
    %           in which the phase diagram will be plotted
    % K0 : [numeric] orbital-orbital coupling for the isotropic 2-channel spin-orbital Kondo(2soK) model
    %           in which the phase diagram will be plotted
    %
    % <Option>
    % 'MeshIntX', .. : [numeric] logarithmic mesh interval for x-axis
    % 'MeshIntY', .. : [numeric] logarithmic mesh interval for y-axis
    %
    % <Output>
    % T-I phase diagram of the isotropic 2soK model for given (J0,K0)

    %% Check inputs
    if ~isnumeric(J0)
        error('ERR: J0 must be a real number')
    end

    if ~isnumeric(K0)
        error('ERR: K0 must be a real number');
    end

    %% Parse options

    while ~isempty(varargin)
        switch varargin{1}
            case 'MeshIntX'
                if ~isnumeric(varargin{2})
                    error('ERR: ''MeshIntX'' must be a number');
                elseif varargin{2} <=0
                    error('ERR: ''MeshIntX'' must be a positive real number');
                else
                    MeshIntX = varargin{2};
                    varargin(1:2) = [];
                end
            case 'MeshIntY'
                if ~isnumeric(varargin{2})
                    error('ERR: ''MeshIntY'' must be a number');
                elseif varargin{2} <=0
                    error('ERR: ''MeshIntY'' must be a positive real number');
                else
                    MeshIntY = varargin{2};
                    varargin(1:2) = [];
                end
        end
    end

    %% Obtain phases and their areas from relevant data
    [I, Phase_range, Phase_name] = PhaseRange_TI(J0,K0,'flatThres',0.1,'platLenThres',0.5);
    [I, Idx] = sort(I,'ascend');
    Phase_range = Phase_range(Idx);
    Phase_name = Phase_name(Idx);

    Phase_area = cell(1,4);     % phase areas on the T-I plane. Each cell element corresponds to Fermi liquid, orbital overscreend, spin overscreened, and fully overscreened phase, respectively
    for it = 1:4
        Phase_area{it} = nan(numel(I), 3);    % phase area on the T-I plane. Each row is [I, minT(I), maxT(I)]
    end
    ExistPhase = false(1,4);    % Variable to check existing phases

    for itD = 1:numel(I)    % for all dataset
        for itP = 1:numel(Phase_name{itD})  % for all phases in each dataset
            switch Phase_name{itD}{itP}
                case 'Fermi Liquid'
                    Phase_area{1}(itD,:) = [I(itD), Phase_range{itD}(itP,:)];
                    ExistPhase(1) = true;
                case 'Orbital Overscreened'
                    Phase_area{2}(itD,:) = [I(itD), Phase_range{itD}(itP,:)];
                    ExistPhase(2) = true;
                case 'Spin Overscreened'
                    Phase_area{3}(itD,:) = [I(itD), Phase_range{itD}(itP,:)];
                    ExistPhase(3) = true;
                case 'Fully Overscreened'
                    Phase_area{4}(itD,:) = [I(itD), Phase_range{itD}(itP,:)];
                    ExistPhase(4) = true;
            end
        end % itP
    end % itD

    %% Create custom colormaps for each phase

    Ngrad = 100;
    Colormap = cell(1,4);
    Colormap(:) = {zeros(Ngrad,3)};
    Colormap{1}(1,:) = [.753,.878,.753];    % Fermi liquid phase
    Colormap{2}(1,:) = [.737,.890,.996];
    Colormap{3}(1,:) = [.961,.757,.675];
    Colormap{4}(1,:) = [.745, .682, .898];
    for itP = 1:4
        for itC = 1:3
            Colormap{itP}(:,itC) = linspace(Colormap{itP}(1,itC),1,Ngrad);
        end
    end

    Cmesh = cell(1,2);
    Cmesh{1} = [-power(10, 1:-1:-14); power(10, -14:1:1)];
    Cmesh{2} = power(10, -16:1:1);

    figure;
    hold on;
    [~,lin2sym_X,sym2lin_X] = SLplot(Cmesh{1},Cmesh{2},'XScale','symlog','YScale','log');
    Cmesh{1} = lin2sym_X(Cmesh{1});
    Cmesh{2} = Cmesh{2};

    %% Define colormap for each phase
    Phase_colormap = cell(1,4);

    for itP = 1:4
        Phase_colormap{itP} = ones(numel(Cmesh{1})-1, numel(Cmesh{2})-1, 3);

        % regions inside the phase boundary
        for it1 = 1:size(Phase_colormap{itP},1)
            for it2 = 1:size(Phase_colormap{itP},2)
                
                Xc = ( Cmesh{1}(it1) + Cmesh{1}(it1+1) ) / 2;
                Yc = sqrt( Cmesh{2}(it2) * Cmesh{2}(it2+1) );

                Inside = false;
                for itD = 1:size(Phase_area{itP},1)-1
                    if ~isnan(Phase_area{itP}(itD,1)) && ~isnan(Phase_area{itP}(itD+1,1))
                        if Xc > lin2sym_X(Phase_area{itP}(itD,1)) && Xc < lin2sym_X(Phase_area{itP}(itD+1,1))
                            if log10(Yc) > Phase_area{itP}(itD,2) && log10(Yc) < Phase_area{itP}(itD,3)
                                Inside = true;
                            end
                        end
                    end
                end % itD
                
                if Inside
                    Phase_colormap{itP}(it1,it2,:) = Colormap{itP}(1,:);
                end

            end % it2
        end % it1

        %{
        % regions outside the phase boundary
        for it1 = 1:size(Phase_colormap{itP},1)
            for it2 = 1:size(Phase_colormap{itP},2)
                
                Xc = ( Cmesh{1}(it1) + Cmesh{1}(it1+1) ) / 2;
                Yc = ( Cmesh{2}(it2) + Cmesh{2}(it1+2) ) / 2;

                Inside = false;
                for itD = 1:numel(Phase_area{itP},1)-1
                    if ~isnan(Phase_area{itP}) && ~isnan(Phase_area{itP+1})
                        if Xc > log10(Phase_area{itP}(itD,1)) && Xc < log10(Phase_area{itP}(itD,1))
                            Inside = true;
                        end
                    end
                end % itD

            end % it2
        end % it1
        %}

    end % itP

    Phase_colormap_comb = ones(numel(Cmesh{1})-1, numel(Cmesh{2})-1, 3);

    for it1 = 1:size(Phase_colormap{itP},1)
            for it2 = 1:size(Phase_colormap{itP},2)
                for itP = 1:4
                    if ~isequal(reshape(Phase_colormap{itP}(it1,it2,:),[1,3]),[1,1,1])
                        Phase_colormap_comb(it1,it2,:) = Phase_colormap{itP}(it1,it2,:);
                        %disp(reshape(Phase_colormap_comb(it1,it2,:),[1,3]));
                        %disp(reshape(Phase_colormap{itP}(it1,it2,:),[1,3]));
                    else
                        %disp('Hi');
                    end
                end
            end
    end

    plotColorgrid(Cmesh,Phase_colormap_comb);
    hold off;

    

end