function [J_O, K_p_O, K_z_O, I_p_O, I_z_O, logD_O] = Solve_RG_Aniso(J0, K_p0, K_z0, I_p0, I_z0, varargin)
    Dstep = -1e-4;
    D0 = 1;
    J = J0;
    K_p = K_p0;
    K_z = K_z0;
    I_p = I_p0;
    I_z = I_z0;
    logD = log(D0)/log(10);

    order = 2;
    
    while ~isempty(varargin)
        switch varargin{1}
            case 'order'
                if isnumeric(varargin{2})
                    order = varargin{2};
                    varargin(1:2) = [];
                else
                    error('ERR: order of perturbation must be a interger greater than 1');
                end
            otherwise
                error(['ERR: unknown option',varargin{1}]);
        end
    end

    function diff = Der(var, J, K_p, K_z, I_p, I_z, order)
        switch order
            case 2
                switch var
                    case 'J'
                        diff = - (J^2 + (2*I_p^2 + I_z^2) / 16 );
                    case 'K_p'
                        diff = - ( K_p*K_z + 3*I_p*I_z / 16 );
                    case 'K_z'
                        diff = - ( K_p^2 + 3*I_p^2 / 16 );
                    case 'I_p'
                        diff = - 1/2 * ( 2*J*I_p + K_p*I_z + K_z*I_p );
                    case 'I_z'
                        diff = - ( J*I_z + K_p*I_p );
                end
            case 3
                switch var
                    case 'J'
                        
                    case 'K_p'
                        
                    case 'K_z'
                        
                    case 'I_p'
                        
                    case 'I_z'
                        
                end
        end
    end

    cnt = 1;
    J_O = [];
    K_p_O = [];
    K_z_O = [];
    I_p_O = [];
    I_z_O = [];
    logD_O = [];

    while(logD > log(1e-16))
        der = zeros(4,5);
        der(1,1) = Der('J', J, K_p, K_z, I_p, I_z, order);
        der(1,2) = Der('K_p', J, K_p, K_z, I_p, I_z, order);
        der(1,3) = Der('K_z', J, K_p, K_z, I_p, I_z, order);
        der(1,4) = Der('I_p', J, K_p, K_z, I_p, I_z, order);
        der(1,5) = Der('I_z', J, K_p, K_z, I_p, I_z, order);
    
        J_tmp = J + (Dstep/2)*der(1,1);
        K_p_tmp = K_p + (Dstep/2)*der(1,2);
        K_z_tmp = K_z + (Dstep/2)*der(1,3);
        I_p_tmp = I_p + (Dstep/2)*der(1,4);
        I_z_tmp = I_z + (Dstep/2)*der(1,5);
        der(2,1) = Der('J', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(2,2) = Der('K_p', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(2,3) = Der('K_z', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(2,4) = Der('I_p', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(2,5) = Der('I_z', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
    
        J_tmp = J + (Dstep/2)*der(2,1);
        K_p_tmp = K_p + (Dstep/2)*der(2,2);
        K_z_tmp = K_z + (Dstep/2)*der(2,3);
        I_p_tmp = I_p + (Dstep/2)*der(2,4);
        I_z_tmp = I_z + (Dstep/2)*der(2,5);
        der(3,1) = Der('J', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(3,2) = Der('K_p', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(3,3) = Der('K_z', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(3,4) = Der('I_p', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(3,5) = Der('I_z', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
    
        J_tmp = J + Dstep*der(3,1);
        K_p_tmp = K_p + Dstep*der(3,2);
        K_z_tmp = K_z + Dstep*der(3,3);
        I_p_tmp = I_p + Dstep*der(3,4);
        I_z_tmp = I_z + Dstep*der(3,5);
        der(4,1) = Der('J', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(4,2) = Der('K_p', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(4,3) = Der('K_z', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(4,4) = Der('I_p', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
        der(4,5) = Der('I_z', J_tmp, K_p_tmp, K_z_tmp, I_p_tmp, I_z_tmp, order);
    
        weights = [1, 2, 2, 1];
        J = J + (Dstep/6)*weights*der(:,1);
        K_p = K_p + (Dstep/6)*weights*der(:,2);
        K_z = K_z + (Dstep/6)*weights*der(:,3);
        I_p = I_p + (Dstep/6)*weights*der(:,4);
        I_z = I_z + (Dstep/6)*weights*der(:,5);
        logD = logD + Dstep;

        if rem(cnt,20) == 0
            J_O = [J_O, J];
            K_p_O = [K_p_O, K_p];
            K_z_O = [K_z_O, K_z];
            I_p_O = [I_p_O, I_p];
            I_z_O = [I_z_O, I_z];
            logD_O = [logD_O, logD];
        end

        if rem(cnt,2e4) == 0
            disp(logD/log(10));
        end
        cnt = cnt + 1;
    end


end