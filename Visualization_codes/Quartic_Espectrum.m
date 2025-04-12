function Espec_Table = Quartic_Espectrum(Lambda_z, Lambda_x, N_iter, Nshow)
    % <Description>
    % Creates a table containing the low-energy spectrum of the 4-electron interaction Kondo model
    % 
    % <Input>
    % Lambda_z : [numeric] Orbital-z interaction strength of the 4-el. int. Kondo model
    % Lambda_x : [numeric] 4-electron interaction strength of the 4-el. int. Kondo model
    % N_iter : [numeric] The NRG iteration step to extract data from (first bath site = 1)
    % N_show : [numeric] The number of lowest energy multiplets to be shown in the table
    %
    % <Output>
    % Espec_Table : [table] Table containing the low-energy spectrum of the 4-electron interaction Kondo model.
    %               The energy levels are normalized so that the ground state energy is 0 and the first excited state energy is 1

    path = [fileparts(mfilename('fullpath')),filesep,'QuarticData_Espectrum',filesep,'K_z=', num2str(2*Lambda_z), '_Q=', num2str(Lambda_x), '_T=1e-24_Nkeep=3000'];
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
    Charge_1 = round(Elev(:,2));
    Charge_2 = round(Elev(:,3));
    Spin_1 = round(Elev(:,4));
    Spin_2 = round(Elev(:,5));

    Espec_Table = table(Energy_levels,Charge_1,Charge_2,Spin_1,Spin_2);
end