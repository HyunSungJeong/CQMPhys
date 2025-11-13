function [J_out,Kp_out,Kz_out,Ip_out,Iz_out,log_D_out] = RG_PoorMan_Aniso(J0,Kp0,Kz0,Ip0,Iz0,varargin)

    h = -1e-4;
    D0 = 1;
    J = J0;
    Kp = Kp0;
    Kz = Kz0;
    Ip = Ip0;
    Iz = Iz0;
    log_D = log(D0)/log(10);
    order = 2;
    dispD = false;
    T = 1e-16;
    
    while ~isempty(varargin)
        switch varargin{1}
            case 'order'
                if isnumeric(varargin{2})
                    order = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: order of perturbation must be a interger greater than 1');
                end

            case 'Temp'
                T = varargin{2};
                varargin(1:2) = [];

            case '-v'
                dispD = true;
                varargin(1) = [];

            otherwise
                error(['ERR: unknown option',varargin{1}]);
        end
    end

    function D = Der(var,J,Kp,Kz,Ip,Iz,order)

        M = 1e4;

        if isinf(J)
            J = sign(J)*M;
        end
        if isinf(Kp)
            Kp = sign(Kp)*M;
        end
        if isinf(Kz)
            Kz = sign(Kz)*M;
        end
        if isinf(Ip)
            Ip = sign(Ip)*M;
        end
        if isinf(Iz)
            Iz = sign(Iz)*M;
        end
        

        switch order
            case 2
                switch var
                    case 'J'
                        D = - (J*J + (1/16)*(2*Ip*Ip + Iz*Iz));
                    case 'Kp'
                        D = - (Kp*Kz + (3/16)*Ip*Iz);
                    case 'Kz'
                        D = - (Kp*Kp + (3/16)*Ip*Ip);
                    case 'Ip'
                        if Ip == 0 && Iz == 0
                            D = 0;
                        else
                            D = - (2*J*Ip + Kp*Iz + Kz*Ip);
                        end
                    case 'Iz'
                        D = - 2*(J*Iz + Kp*Ip);
                end
            case 3
                switch var
                    case 'J'

                    case 'Kp'

                    case 'Kz'

                    case 'Ip'
                        
                    case 'Iz'
                end
        end
    end

    cnt = 1;
    J_out = [];
    Kp_out = [];
    Kz_out = [];
    Ip_out = [];
    Iz_out = [];
    log_D_out = [];

    while(log_D > log(T))
        der = zeros(4,5);
        der(1,1) = Der('J',J,Kp,Kz,Ip,Iz,order);
        der(1,2) = Der('Kp',J,Kp,Kz,Ip,Iz,order);
        der(1,3) = Der('Kz',J,Kp,Kz,Ip,Iz,order);
        der(1,4) = Der('Ip',J,Kp,Kz,Ip,Iz,order);
        der(1,5) = Der('Iz',J,Kp,Kz,Ip,Iz,order);
    
        J_temp = J + (h/2)*der(1,1);
        Kp_temp = Kp + (h/2)*der(1,2);
        Kz_temp = Kz + (h/2)*der(1,3);
        Ip_temp = Ip + (h/2)*der(1,4);
        Iz_temp = Iz + (h/2)*der(1,5);

        der(2,1) = Der('J',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(2,2) = Der('Kp',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(2,3) = Der('Kz',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(2,4) = Der('Ip',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(2,5) = Der('Iz',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
    
        J_temp = J + (h/2)*der(2,1);
        Kp_temp = Kp + (h/2)*der(2,2);
        Kz_temp = Kz + (h/2)*der(2,3);
        Ip_temp = Ip + (h/2)*der(2,4);
        Iz_temp = Iz + (h/2)*der(2,5);

        der(3,1) = Der('J',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(3,2) = Der('Kp',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(3,3) = Der('Kz',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(3,4) = Der('Ip',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(3,5) = Der('Iz',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
    
        J_temp = J + h*der(3,1);
        Kp_temp = Kp + h*der(3,2);
        Kz_temp = Kz + h*der(3,3);
        Ip_temp = Ip + h*der(3,4);
        Iz_temp = Iz + h*der(3,5);

        der(4,1) = Der('J',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(4,2) = Der('Kp',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(4,3) = Der('Kz',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(4,4) = Der('Ip',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
        der(4,5) = Der('Iz',J_temp,Kp_temp,Kz_temp,Ip_temp,Iz_temp,order);
    
        weights = [1, 2, 2, 1];
        J = J + (h/6)*weights*der(:,1);
        Kp = Kp + (h/6)*weights*der(:,2);
        Kz = Kz + (h/6)*weights*der(:,3);
        Ip = Ip + (h/6)*weights*der(:,4);
        Iz = Iz + (h/6)*weights*der(:,5);
        
        log_D = log_D + h;

        N = 1e4;
        if isinf(J)
            J = sign(J)*N;
        end
        if isinf(Kp)
            Kp = sign(Kp)*N;
        end
        if isinf(Kz)
            Kz = sign(Kz)*N;
        end
        if isinf(Ip)
            Ip = sign(Ip)*N;
        end
        if isinf(Iz)
            Iz = sign(Iz)*N;
        end

        if rem(cnt,20) == 0
            J_out = [J_out, J];
            Kp_out = [Kp_out, Kp];
            Kz_out = [Kz_out, Kz];
            Ip_out = [Ip_out, Ip];
            Iz_out = [Iz_out, Iz];
            log_D_out = [log_D_out, log_D];
        end

        if rem(cnt,2e4) == 0 && dispD
            disp(log_D/log(10));
        end
        cnt = cnt + 1;
    end

end