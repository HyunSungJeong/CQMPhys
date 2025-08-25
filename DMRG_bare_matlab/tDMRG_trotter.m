function [t_steps, MPS_final, O_vals, dw] = tDMRG_trotter(MPS, Op_loc, H_NN, Nkeep, dt, Tmax)
    % <Description>
    % Performs tDMRG to the given MPS according to the given NN Hamiltonian
    %
    % <Input>
    % MPS : [cell] MPS tensors of the initial state. MPS{n} is the n-th MPS tensor
    %
    % Op_loc : [cell] local operators to compute the expectation value at
    %                 each time step of tDMRG
    %
    % H_NN : [cell] nearest-neighbor interations in the system Hamiltonian.
    %               H_NN{n} is the rank-4 interaction corresponding to the
    %               interaction Hamiltonian between n-th and (n+1)-th site
    %
    % Nkeep : [numeric] maximum bond dimension of the MPS to be kept
    %
    % dt : [numeric] time step size to be used in the trotterization
    %
    % Tmax : [numeric] maximum time
    %
    % <Options>
    % 
    %
    % <Output>
    % t_steps : [numeric vector] row vector of tDMRG time steps
    %
    % MPS_final : [cell] MPS at time Tmax
    %
    % O_vals : [cell] Local operator expectation values at each time step.
    %                 O_vals{itN}(itT) is the expectation value Op_loc{itN} at time t_steps(itT)
    %
    % dw : [numeric] matrix of discarded weights.
    %                dw(m,n) is the discarded weight at m-th trotterization step, between n-th and (n+1)-th site

    
end