function [J_out,K_out,I_out,log_D_out] = RG_PoorMan(J0,K0,I0,varargin)

    h = -1e-4;
    D0 = 1;
    J = J0;
    K = K0;
    I = I0;
    log_D = log(D0)/log(10);
    order = 2;
    Rho = 1;
    dispD = false;
    
    while ~isempty(varargin)
        switch varargin{1}
            case 'order'
                if isnumeric(varargin{2})
                    order = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: order of perturbation must be a interger greater than 1');
                end

            case 'Rho'
                if isnumeric(varargin{2})
                    Rho = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: ''Hyb'' must be a positive real number');
                end

            case '-v'
                dispD = true;
                varargin(1) = [];

            otherwise
                error(['ERR: unknown option',varargin{1}]);
        end
    end

    function D = Der(var,J,K,I,Rho,order)
        switch order
            case 2
                switch var
                    case 'J'
                        D = - (J*J + (3/16)*I*I);
                    case 'K'
                        D = - (K*K + (3/16)*I*I);
                    case 'I'
                        D = - 2*(J+K)*I;
                end
            case 3
                switch var
                    case 'J'
                        D = - (1 - J)*(J*J + (3/16)*I*I);
                    case 'K'
                        D = - (1 - K)*(K*K + (3/16)*I*I);
                    case 'I'
                        D = - ((2*J + 2*K - J*J - K*K)*I - (1/8)*I*I*I);
                end
        end

        D = Rho*D;
    end

    cnt = 1;
    J_out = [];
    K_out = [];
    I_out = [];
    log_D_out = [];

    while(log_D > log(1e-24))
        der = zeros(4,3);
        der(1,1) = Der('J',J,K,I,Rho,order);
        der(1,2) = Der('K',J,K,I,Rho,order);
        der(1,3) = Der('I',J,K,I,Rho,order);
    
        J_temp = J + (h/2)*der(1,1);
        K_temp = K + (h/2)*der(1,2);
        I_temp = I + (h/2)*der(1,3);
        der(2,1) = Der('J',J_temp,K_temp,I_temp,Rho,order);
        der(2,2) = Der('K',J_temp,K_temp,I_temp,Rho,order);
        der(2,3) = Der('I',J_temp,K_temp,I_temp,Rho,order);
    
        J_temp = J + (h/2)*der(2,1);
        K_temp = K + (h/2)*der(2,2);
        I_temp = I + (h/2)*der(2,3);
        der(3,1) = Der('J',J_temp,K_temp,I_temp,Rho,order);
        der(3,2) = Der('K',J_temp,K_temp,I_temp,Rho,order);
        der(3,3) = Der('I',J_temp,K_temp,I_temp,Rho,order);
    
        J_temp = J + h*der(3,1);
        K_temp = K + h*der(3,2);
        I_temp = I + h*der(3,3);
        der(4,1) = Der('J',J_temp,K_temp,I_temp,Rho,order);
        der(4,2) = Der('K',J_temp,K_temp,I_temp,Rho,order);
        der(4,3) = Der('I',J_temp,K_temp,I_temp,Rho,order);
    
        weights = [1, 2, 2, 1];
        J = J + (h/6)*weights*der(:,1);
        K = K + (h/6)*weights*der(:,2);
        I = I + (h/6)*weights*der(:,3);
        log_D = log_D + h;

        if rem(cnt,20) == 0
            J_out = [J_out, J];
            K_out = [K_out, K];
            I_out = [I_out, I];
            log_D_out = [log_D_out, log_D];
        end

        if rem(cnt,2e4) == 0 && dispD
            disp(log_D/log(10));
        end
        cnt = cnt + 1;
    end

end