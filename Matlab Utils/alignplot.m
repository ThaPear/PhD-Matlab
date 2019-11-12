function alignplot(hFig, Nx, Ny, Px, Py, screen)
    if(verLessThan('matlab','9.3'))
        % This function does not work in 2016b or older, so don't try.
        return;
    end
    if(nargin < 6)
        screen = 2;
    end
    if(nargin < 5 || isempty(Py))
        Py = floor((Px-1)/Nx)+1;
        Px = mod(Px-1, Nx)+1;
    end
    
    Px = mod(Px-1, Nx)+1;
    Py = mod(Py-1, Ny)+1;
    
    gr = groot;
    x0 = gr.MonitorPositions(screen, 1);
    y0 = gr.MonitorPositions(screen, 2);
    szx = gr.MonitorPositions(screen, 3);
    szy = gr.MonitorPositions(screen, 4);
    
    % Subtract the start bar.
%     if(screen == 1)
        szx = szx - 70;
%     end
    
    hFig.OuterPosition = [x0+(Px-1)/Nx*szx, y0+(Ny-(Py))/Ny*szy, ...
                          szx/Nx, szy/Ny];
                      
    
end
