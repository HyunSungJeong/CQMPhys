function test(parfn)
    partot = job_func_preamble(parfn);

    [PE] = loadvar(partot, ...
        {'PE'}, ...
        {[]});
    text = testfunc();
    save('/home/hyunsung/test.mat','text');
end