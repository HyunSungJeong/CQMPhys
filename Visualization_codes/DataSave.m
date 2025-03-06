clear;


for it = (11:14)
    Cpl_Const = [0.1,0.3,power(10,-it)];
    key = key_enc(Cpl_Const);
    Info.num_phase = 3;
    Info.phase = {'Unscreened','NFL1','NFL2'};
    Info.boundaries = [-3.015, 15.76];
    Phase_Info = containers.Map(key, 1);
end

function key1_char =  key_enc(key1)
    [m,n] = size(key1);
    key1_char=repmat({char(ones(1,16*n))},m,1);
    tmp = num2hex(key1);
    for ii=1:m  
    key1_char{ii,:} = reshape((tmp(ii:m:end,:).'),1,16*n);
    end
end
