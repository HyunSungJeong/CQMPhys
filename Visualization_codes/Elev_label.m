function Elev_labeled = Elev_label(Etot, Qtot, N, N_elev)
    % <Description>
    % <Input>
    % Etot: energy level data from the output of 'plotE' function
    % Qtot: quantum number data from the output of 'plotE' function
    % N : the iteration step to analyze (counting from 0)
    % N_elev : the number of lowest energy levels to be shown with their quantum numbers
    %
    % <Output>
    % Elev_labeled : [1x2 cell array]
    %               the lowest N_elev energy levels with their quantum numbers from NRG iteration #N
    %               the 1st column is the lowest energy levels in increasing order,
    %               and the 2nd-5th columns are the corresponding quantum numbers for each energy levels

    N = N+1;
    
    if ~ismember(N, 1:numel(Etot))
        error('ERR: N must be a nonnegative integer that does not exceed Wilson chain length');
    end

    Elevs = [];
    for it = 1:numel(Etot{N})
        Elevs = cat(1, Elevs, sort(Etot{N}{it}(:)));
    end

    Qnums = [];
    for it1 = 1:size(Qtot{N}{1},1)
        for it2 = 1:numel(Etot{N}{it1})
            Qnums = cat(1, Qnums, [Qtot{N}{1}(it1,:), it2] );
        end
    end

    [Elevs, idx] = sort(Elevs);
    Qnums = Qnums(idx,:);

    Elev_labeled = [Elevs(1:N_elev), Qnums(1:N_elev,:)];
end