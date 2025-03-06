function G_o = Gfk_o(epsilon, omega, SE_o, U, V, t_0, t, t_p, muc, ieta, idx1, idx2)    % lattice GF. in k-space, odd sector
    % < Input >
    % epsilon : [numeric] noninteracting band energy
    % omega : [numeric] frequency
    % SE_o : [numeric] 2x2 matrix valued self-energy in the even sector, at frequency omega
    % U : [numeric] on-site Coulomb repulsion
    % V : [numeric] inter-site hopping amplitude
    % t_0 : [numeric] intra-cell hopping amplitude t_0
    % t : [numeric] intra-cell hopping amplitude t
    % t_p : [numeric] intra-cell hopping amplitude t'
    % muc : [numeric] chemical potential
    % ieta : [numeric] convergence generating factor
    % idx1 : [numeric] index of the GF. along its first dimension (1<=idx1<=2)
    % idx2 : [numeric] index of the GF. along its second dimension (1<=idx2<=2)

    % < Output >
    % G_o : [numeric] (idx1, idx2) component of the lattice GF. of the odd sector in k-space

    G_o = zeros(1,numel(epsilon));
    omega_p = omega + ieta;     % omega with infinitesimal imaginary part

    for it = 1:numel(epsilon)
        G = [omega_p + muc + U/2 - epsilon(it), -(t-t_p);
            -(t-t_p), omega_p + muc + U/2];

        G = G - SE_o;
        G = inv(G);
        G_o(it) = G(idx1,idx2);
    end
end