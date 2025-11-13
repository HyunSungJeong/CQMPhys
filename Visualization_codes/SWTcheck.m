function [I, J0, K_perp, K_z, mu] = SWTcheck(U, J, JL, Hyb)
    
    % Compute candidate bounds for mu
    expr1 = U - 2*J + JL/4;
    expr2 = U - (3/4)*JL;
    expr3 = U + JL/4;
    
    minExpr = min(min(expr1, expr2), expr3);
    
    % Define mu from condition in the problem
    mu = 0.5 * (minExpr - JL/4);
    gap = 0.5 * (minExpr + JL/4);
    
    % Check inequality alpha:  -JL/4 < mu < minExpr
    cond1 = (-JL/4 < mu);
    cond2 = (mu < minExpr);
    valid = cond1 & cond2;
    
    % Define I
    I = 1 ./ (U - 2*J + JL/4 - mu) + ...
    1 ./ (U - 3*JL/4 - mu) + ...
    2 ./ (JL/4 + mu);

    J0 = -1 ./ (U-2*J + JL/4 - mu) ...
        +2 ./ (U - 3*JL/4 - mu) ...
        +1 ./ (U + JL/4 - mu) ...
        +2 ./ (JL/4 + mu);
    J0 = J0 ./ 2;
    
    K_perp = 3 ./ (U-2*J + JL/4 - mu) ...
             -1 ./ (U + JL/4 - mu) ...
             +2 ./ (JL/4 + mu);
    K_perp = K_perp ./ 2;
    
    K_z = 3 ./ (U-2*J + JL/4 - mu) ...
          -2 ./ (U - 3*JL/4 - mu) ...
          +1 ./ (U + JL/4 - mu) ...
          +2 ./ (JL/4 + mu);
    K_z = K_z ./ 2;
    
    I = Hyb*I;  J0 = Hyb*J0;    K_perp = Hyb*K_perp;   K_z = Hyb*K_z;

    if ~valid
        disp('WRN: The 2oAH model cannot be projected to the 2soK model');
        I = NaN; J0 = NaN; K_perp = NaN; K_z = NaN; mu = NaN;
    end
end