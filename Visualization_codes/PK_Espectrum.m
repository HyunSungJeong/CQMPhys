function Espec = PK_Espectrum(Lambda_z, Lambda_x, N_iter, Nshow)
    % <Description>
    % Creates a table containing the low-energy spectrum of the pair-Kondo model
    % 
    % <Input>
    % Lambda_z : [numeric] Orbital-z interaction strength of the pair-Kondo model
    % Lambda_x : [numeric] 4-electron interaction strength of the pair-Kondo model
    % N_iter : [numeric] The NRG iteration step to extract data from (counting from 0)
    % N_show : [numeric] The number of lowest energy multiplets to be shown in the table
    %
    % <Output>
    % Espec : [table] Table containing the low-energy spectrum of the pair-Kondo model.
    %               The energy levels are normalized so that the ground state energy is 0 and the first excited state energy is 1

    path = [fileparts(mfilename('fullpath')),filesep,'PKData_Espectrum',filesep,'lambda_z=', num2str(Lambda_z), '_lambda_x=', num2str(Lambda_x), '_T=1e-24_Nkeep=3000'];
    Etot = load([path, '\Etot.mat']);
    Etot = Etot.Etot;

    Qtot = load([path, '\Qtot.mat']);
    Qtot = Qtot.Qtot;

    Elev = Elev_label(Etot,Qtot,N_iter,Nshow);
    
    it = 1;
    while Elev(it,1) < 0.01
        it = it + 1;
        E1 = Elev(it,1);
    end

    Energy_levels = Elev(:,1)/E1;
    U1charge_1 = round(Elev(:,2));
    U1charge_2 = round(Elev(:,3));
    SU2spin_1 = round(Elev(:,4))/2;
    SU2spin_2 = round(Elev(:,5)).2;

    Espec = table(Energy_levels,U1charge_1,U1charge_2,SU2spin_1,SU2spin_2);
end