function resultdir = resultdir_cst()
    resultdir = 'e:/ Simulations/Matlab-Generated';
    
    [~, hostname] = system('hostname'); hostname = strsplit(hostname, '\n');
    switch(lower(hostname{1}))
        case 'srv539'
            path = 'e:/data/Sander/Matlab-Generated';
        case 'tud211735'
            path = 'e:/ Simulations/Matlab-Generated';
        otherwise
            resultdir = resultdir_matlab;
    end
end