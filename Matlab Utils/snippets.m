%% Equal axis scales.
    set(gca,'dataAspectRatio',[1 1 1]);
    
%% Progress bar.
progressbar = waitbar(0, 'Calculating...',...
                'Name','Progress');
for(i = 1:100)
        waitbar(i/100, progressbar, ...
            ['Calculating... ', num2str(f/max(fs)*100, '%.1f'), '%']);
end
close(progressbar);