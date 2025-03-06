function G_e = Gfk_e(epsilon, omega, SE_e, U, V, t_0, t, t_p, muc, ieta, idx1, idx2)    % lattice GF. in k-space, even sector
    % < Input >
    % epsilon : [numeric] noninteracting band energy
    % omega : [numeric] frequency
    % SE_e : [numeric] 3x3 matrix valued self-energy in the even sector, at frequency omega
    % U : [numeric] on-site Coulomb repulsion
    % V : [numeric] inter-site hopping amplitude
    % t_0 : [numeric] intra-cell hopping amplitude t_0
    % t : [numeric] intra-cell hopping amplitude t
    % t_p : [numeric] intra-cell hopping amplitude t'
    % muc : [numeric] chemical potential
    % ieta : [numeric] convergence generating factor
    % idx1 : [numeric] index of the GF. along its first dimension (1<=idx1<=3)
    % idx2 : [numeric] index of the GF. along its second dimension (1<=idx2<=3)

    % < Output >
    % G_e : [numeric] (idx1, idx2) component of the lattice GF. of the even sector in k-space

    G_e = zeros(1,numel(epsilon));
    omega_p = omega + ieta;     % omega with infinitesimal imaginary part

    for it = 1:numel(epsilon)
        G = [omega_p + muc + U/2 - epsilon(it), -(t+t_p), 0;
            -(t+t_p), omega_p + muc + U/2, -sqrt(2)*t_0;
            0, -sqrt(2)*t_0, omega_p + muc + U/2];

        G = G - SE_e;
        G = inv(G);
        G_e(it) = G(idx1,idx2);
    end
end