
classdef NPort < handle
    properties
        ports
        tline
        
        STE
        STM
        S_fs
        S_th
        S_ph
    end
    properties(Dependent)
        N
    end
    methods
        function this = NPort(port, terminatedtline)
            % Connects the given port to the input terminals of the terminated tline.
            if(isa(terminatedtline, 'TerminatedTLine'))
                this.tline = TLine({port, terminatedtline.tline, terminatedtline.terminator});
            elseif(isa(terminatedtline, 'TLine'))
                this.tline = TLine({port, terminatedtline});
            else
                error('Unknown type ''%s''.\n', class(terminatedtline));
            end
            
            % Find all ports in the network.
            this.ports = NPort.FindPorts(this.tline);
            
            % Sort the ports by their number.
            portnumbers = [this.ports.number];
            [~, I] = sort(portnumbers);
            this.ports = this.ports(I);
            
%             dispex('Created NPort with %i ports.\n', this.N);
        end
        function n = get.N(this)
            n = length(this.ports);
        end
        function InitializeSMatrix(this, fs, th, ph)
            if(th == 0)
                error('Cannot have theta 0.');
            end
            if(~isempty(this.S_fs) && (th ~= this.S_th || ph ~= this.S_ph))
                this.STE = [];
                this.STM = [];
                this.S_fs = [];
                this.S_th = [];
                this.S_ph = [];
            end
            newfs = setdiff(fs, this.S_fs);
            if(length(newfs) < 1)
                return;
            end
            dispex('S Matrix: Calculating for %i frequencies, %i ports.\n', length(fs), this.N);
            tc = tic;
            
            Np = this.N;
            Nf = length(newfs);
            
            %% Initialize progress bar
            hDataQueue = parallel.pool.DataQueue;
            hWaitbar = waitbar(0, {'0% Initializing...', ''});
            afterEach(hDataQueue, @updateWaitbar);
            progress = -1; send(hDataQueue, nan);
            
            this.S_fs = [this.S_fs, newfs];
            
            tline_ = this.tline;
            
            newSTE = zeros(Np, Np, Nf);
            newSTM = zeros(Np, Np, Nf);
            parfor(fi = 1:Nf) % parfor
                f = newfs(fi);
                
                newSTEfi = zeros(Np, Np);
                newSTMfi = zeros(Np, Np);
                for(i = 1:Np)
                    port1 = this.ports(i);
                    for(j = i+1:Np)
                        port2 = this.ports(j);
                        % Reorganize the TLine such that it represents the TLine from port1 to port2.
                        tlineij = NPort.ReorganizeTLine(tline_, port1, port2);
                        
                        % Calculate the S-matrix of this network.
                        
                        [k0, kx0, ky0, ~] = k(f, 1, th, ph);
                        kr = sqrt(kx0.^2 + ky0.^2);
                        % ABCD Matrix.
                        ABCDTE = tlineij.GetABCD(1, f, k0, kr);
                        ABCDTM = tlineij.GetABCD(0, f, k0, kr);
                        
                        STE_ = ABCD2S(ABCDTE, port1.impedance, port2.impedance);
                        STM_ = ABCD2S(ABCDTM, port1.impedance, port2.impedance);
                        
                        if(newSTEfi(i, i) ~= 0)
                            if(abs(newSTEfi(i, i) - STE_.s11) > 1e-12)
                                % Sii for i>1 has already been calculated through Sjj, so check if
                                % the result is same.
                                error('Found different S%i%i, error is %g.', i, i, abs(newSTEfi(i, i) - STE_.s11));
                            end
                        end
                        % Since we set the other ports to open, the Z matrix for these ports can
                        % simply be used for the big matrix.
                        newSTEfi(i, i) = STE_.s11;                  newSTMfi(i, i) = STM_.s11;
                        newSTEfi(i, j) = STE_.s12;                  newSTMfi(i, j) = STM_.s12;
                        newSTEfi(j, i) = STE_.s21;                  newSTMfi(j, i) = STM_.s21;
                        newSTEfi(j, j) = STE_.s22;                  newSTMfi(j, j) = STM_.s22;
                    end
                end
                % Store the resulting matrices.
                newSTE(:, :, fi) = newSTEfi;
                newSTM(:, :, fi) = newSTMfi;
                % Update progress bar.
                send(hDataQueue, nan);
            end
            % Store values in this object.
            for(fi = 1:Nf)
                f = newfs(fi);
                fii = find(this.S_fs == f);
                this.STE(:, :, fii) = newSTE(:, :, fi);
                this.STM(:, :, fii) = newSTM(:, :, fi);
            end
            this.S_th = th;
            this.S_ph = ph;
            
            function updateWaitbar(~)
                progress = progress + 1;
                waitbar(progress/Nf, hWaitbar, ...
                    {sprintf('%.1f%% Initializing S Matrix...', progress/Nf/2*100), ...
                     sprintf('%i/%i frequencies done.', progress, Nf)});
            end
            delete(hWaitbar);
            
            dt = toc(tc);
            dispex('S Matrix: Completed in %s.\n', fancyduration(dt));
        end
        function [ZTE, ZTM, portnumbers] = GetZMatrix(this, fs, th, ph)
            if(th == 0)
                error('Cannot have theta 0.');
            end
            error('Not implemented.');
            % Ensure Z matrix is initialized for the given frequencies and angles.
            this.InitializeSMatrix(fs, th, ph);
            
            Np = this.N;
            Nf = length(fs);
            
            ZTE = zeros(Np, Np, Nf);
            ZTM = zeros(Np, Np, Nf);
            for(fi = 1:Nf)
                f = fs(fi);
                fii = find(this.S_fs == f);
                
                
                
            end
            % The sequence of the ports in the Z-matrix.
            portnumbers = [this.ports.number];
        end
        function [zref, portnumbers] = GetZRef(this)
            zref = [this.ports.impedance];
            % The sequence of the ports in the Z-vector.
            portnumbers = [this.ports.number];
        end
        function [Zref, portnumbers] = GetZRefMatrix(this)
            Zref = diag([this.ports.impedance]);
            % The sequence of the ports in the Z-matrix.
            portnumbers = [this.ports.number];
        end
        function [Gref, portnumbers] = GetGRefMatrix(this)
            zref = this.GetZRef();
            Gref = diag(1./sqrt(abs(zref)));
            % The sequence of the ports in the G-matrix.
            portnumbers = [this.ports.number];
        end
        function [YTE, YTM, portnumbers] = GetYMatrix(this, fs, th, ph)
            if(th == 0)
                error('Cannot have theta 0.');
            end
            error('Not implemented.');
            % Ensure Z matrix is initialized for the given frequencies and angles.
            this.InitializeSMatrix(fs, th, ph);
            
            Np = this.N;
            Nf = length(fs);
            
            YTE = zeros(Np, Np, Nf);
            YTM = zeros(Np, Np, Nf);
            for(fi = 1:Nf)
                f = fs(fi);
                fii = find(this.S_fs == f);
                
                
                
            end
        end
        function [STE_, STM_, portnumbers] = GetSMatrix(this, fs, th, ph)
            if(th == 0)
                error('Cannot have theta 0.');
            end
            % Ensure Z matrix is initialized for the given frequencies and angles.
            this.InitializeSMatrix(fs, th, ph);
            
            Np = this.N;
            Nf = length(fs);
            
            STE_ = zeros(Np, Np, Nf);
            STM_ = zeros(Np, Np, Nf);
            for(fi = 1:Nf)
                f = fs(fi);
                % Ensure the frequency ordering matches the given frequency list.
                fii = find(this.S_fs == f);
                STE_(:, :, fi) = this.STE(:, :, fii);
                STM_(:, :, fi) = this.STM(:, :, fii);
            end
            % The sequence of the ports in the S-matrix.
            portnumbers = [this.ports.number];
        end
        function [ZinTE, ZinTM, portnumbers] = GetInputImpedance(this, fs, th, ph, excitation)
            if(th == 0)
                error('Cannot have theta 0.');
            end
            % Ensure Z matrix is initialized for the given frequencies and angles.
            this.InitializeSMatrix(fs, th, ph);
            
            dispex('Active Z: Calculating for %i frequencies, %i ports.\n', length(fs), this.N);
            tc = tic;
            
            if(~iscolumn(excitation))
                warning('Expected column excitation, transposing.');
                excitation = excitation.';
            end
            
            Np = this.N;
            Nf = length(fs);

            ZinTE = zeros(Np, Nf);
            ZinTM = zeros(Np, Nf);
            for(fi = 1:Nf)
                f = fs(fi);
                
                % Retrieve the Z matrix for this frequency
                fii = find(this.S_fs == f);
                ZTE_ = this.STE(:, :, fii); ZTM_ = this.STM(:, :, fii);
                % Reference loads.
                Zref = this.GetZRefMatrix();
                % Add the loads to the mutual impedances.
                ZpTE = ZTE_ + Zref;         ZpTM = ZTM_ + Zref;
                % Determine current through the elements.
                iTE = ZpTE\excitation;      iTM = ZpTM\excitation;
                % Determine voltage on the elements themselves.
                vTE = ZTE_ * iTE;           vTM = ZTM_ * iTM;

                ZinTE(:, fi) = vTE./iTE;    ZinTM(:, fi) = vTM./iTM;
            end

            dt = toc(tc);
            dispex('Active Z: Completed in %s.\n', fancyduration(dt));
            % The sequence of the ports in the S-matrix.
            portnumbers = [this.ports.number];
        end
    end
    methods(Static)
        function ports = FindPorts(tline)
            ports = {};
            for(eli = 1:length(tline.elements))
                el = tline.elements{eli};
                if(isa(el, 'Port'))
                    ports = [ports, el];
                elseif(isprop(el, 'elements'))
                    ports = [ports, NPort.FindPorts(el)];
                elseif(isprop(el, 'impedance'))
                    if(isa(el.impedance, 'Port')) % Shunt / Series Port
                        ports = [ports, el.impedance];
                    elseif(isprop(el.impedance, 'elements')) % Shunt / Series Tline
                        ports = [ports, NPort.FindPorts(el.impedance)];
                    end
                end
            end
        end
        function tline = SetPortsToOpen(tline)
            % Sets the impedance of all Port objects in the given tline to Inf.
            for(eli = 1:length(tline.elements))
                el = tline.elements{eli};
                if(isa(el, 'Port'))
                    tline.elements{eli}.impedance = Inf;
                elseif(isprop(el, 'elements'))
                    tline.elements{eli} = NPort.SetPortsToOpen(el);
                elseif(isprop(el, 'impedance'))
                    if(isa(el.impedance, 'Port')) % Shunt / Series Port
                        tline.elements{eli}.impedance.impedance = Inf;
                    elseif(isprop(el.impedance, 'elements')) % Shunt / Series Tline
                        el.impedance = NPort.SetPortsToOpen(el.impedance);
                        tline.elements{eli} = el;
                    end
                end
            end
        end
        function tline = ReorganizeTLine(tline, port1, port2)
            % Reorganizes the TLine such that port1 can be connected to the start of it and to the
            % end of it.
            port1found = 0;
            port2found = 0;
            
            
%             if(port2 == baseport)
%                 % If port2 is the base port, we need to flip the TLine and add port2 to the end.
%                 tline = tline.Flip();
%                 tline.elements{end+1} = port2;
%                 if(tline.elements{1} == port1)
%                     % port1 is already in front, so remove it and return.
%                     tline.elements(1) = [];
%                     return;
%                 end
%                 
%                 % port2 is found
%                 port2found = 1;
%             end
%             if(port1 == baseport)
%                 % We only have to ensure port2 is at the end of the TLine.
%                 if(tline.elements{end} == port2)
%                     % port2 is already at the end, so we're done
%                     tline.elements(end) = [];
%                     return;
%                 end
%                 port1found = 1;
%             end
            while(~port1found || ~port2found)
                % Flatten the t-line.
                tline = tline.Flatten();
                for(eli = 1:length(tline.elements))
                    el = tline.elements{eli};
                    if(isa(el, 'Port'))
                        if(~port1found && el == port1)
                            if(eli == 1)
%                                 dispex('1.\n');
                                % Remove the port from the line.
                                tline.elements(1) = [];
                                % Found port 1.
                                port1found = 1;
                                break;
                            else
%                                 dispex('2.\n');
                                % Move the part of the tline before port1 into a shunt.
                                newshunt = tline.elements(1:eli-1);     % The part to move to the shunt.
                                newshunt = newshunt.Flip();             % Flip so that the end of the line connects at the shunt.
                                tline.elements(1:eli) = [];             % Remove the part of the tline that is to go into the shunt.
                                                                        % The port is also removed.
                                tline.elements = [{Shunt(TLine(newshunt))}, tline.elements]; % Recombine the tline.
                                % Found port 1.
                                port1found = 1;
                                break;
                            end
                        elseif(~port2found && el == port2)
                            if(eli == length(tline.elements))
%                                 dispex('3.\n');
                                % Remove the port from the line.
                                tline.elements(end) = [];
                                % If we haven't found port1 yet, it's not in the tline.
                                if(~port1found)
                                    error('Could not find port 1.');
                                end
                                % Found port 2.
                                port2found = 1;
                                break;
                            else
%                                 dispex('4.\n');
                                % If we're not at the end of the line, move the remaining part of the
                                % tline into a shunt.
                                newshunt = tline.elements(eli+1:end);   % The part to move to the shunt.
                                tline.elements(eli:end) = [];           % Remove the part of the tline that is to go into the shunt.
                                                                        % The port is also removed.
                                tline.elements = [tline.elements, {Shunt(TLine(newshunt))}]; % Recombine the tline.
                                % Found port 2.
                                port2found = 1;
                                break;
                            end
                        end
                    elseif(isa(el, 'Shunt'))
                        if(~port1found && NPort.ContainsPort(el.impedance, port1))
                            % If the impedance of the shunt contains the port it must be a TLine.
                            % Swap the part of the TLine we've handled with the TLine in the Shunt.
                            newshunt = tline.elements(1:eli-1);     % The part to move to the shunt.
                            tline.elements(1:eli) = [];             % Remove it from the tline.
                            oldshunt = el.impedance;                % The part that was in the shunt.
                            oldshunt = oldshunt.Flip();             % Flip so that the start of the shunt connects to the line.
% Commented because it should already be flattened.
%                             if(isprop(oldshunt, 'elements'))
%                                 oldshunt = oldshunt.Flatten();      % Flatten it.
%                             end
                            el.impedance = TLine(newshunt).Flip();  % Move the newshunt into the shunt.
                            if(~isempty(newshunt))
                                tline.elements = [{oldshunt}, {el}, tline.elements]; % Recombine the tline.
                            else
                                tline.elements = [{oldshunt}, tline.elements]; % Recombine the tline.
                            end
                            % Restart iteration.
                            break;
                        end
                        if(~port2found && NPort.ContainsPort(el.impedance, port2))
                            % If the impedance of the shunt contains the port it must be a TLine.
                            % Swap the rest of the TLine with the Tline in the Shunt.
                            newshunt = tline.elements(eli+1:end);   % The part to move to the shunt.
                            tline.elements(eli:end) = [];           % Remove it from the tline.
                            oldshunt = el.impedance;                % The part that was in the shunt.
                            el.impedance = TLine(newshunt);         % Move the newshunt into the shunt.
                            if(~isempty(newshunt))
                                tline.elements = [tline.elements, {el}, {oldshunt}]; % Recombine the tline.
                            else
                                tline.elements = [tline.elements, {oldshunt}]; % Recombine the tline.
                            end
                            % Restart iteration.
                            break;
                        end
                    elseif(isa(el, 'Series'))
                        if(~port1found && NPort.ContainsPort(el.impedance, port1))
                            error('Series not implemented.');
                        end
                        if(~port2found && NPort.ContainsPort(el.impedance, port2))
                            error('Series not implemented.');
                        end
                    elseif(strcmp(class(el), 'TLine')) %#ok<STISA> 'isa' also accepts subclasses, which we don't want.
                        error('TLine should already be flattened.');
                    elseif(strcmp(class(el), 'TerminatedTLine')) %#ok<STISA> 'isa' also accepts subclasses, which we don't want.
                        error('TerminatedTLine should be flattened into TLine.');
                    end
                    if(eli == length(tline.elements))
                        % We've reached the end of the line. Check if we've found both ports.
                        if(~port1found && ~port2found)
                            error('Could not find either port %i or %i.', port1.number, port2.number);
                        elseif(~port1found)
                            error('Could not find port %i.', port1.number);
                        elseif(~port2found)
                            error('Could not find port %i.', port2.number);
                        else
                            dispex('Port %i and %i found.\n', port1.number, port2.number);
                        end
                    end
                end
            end
%             dispex('Port %i and %i found.\n', port1.number, port2.number);
        end
        function boolean = ContainsPort(tline, port)
            
            boolean = 0;
            
            if(isa(tline, 'Port'))
                boolean = (tline == port);
                return;
            end
            if(~isa(tline, 'TLine'))
                return;
            end
            for(eli = 1:length(tline.elements))
                el = tline.elements{eli};
                if(isa(el, 'Port'))
                    % Check if it's the port we're looking for.
                    if(el == port)
                        boolean = 1;
                        return;
                    end
                elseif(isa(el, 'Shunt') || isa(el, 'Series'))
                    % Check if the shunt/series impedance contains the port we're looking for.
                    if(NPort.ContainsPort(el.impedance, port))
                        boolean = 1;
                        return;
                    end
                elseif(isa(el, 'TLine'))
                    % Check if the transmission line contains the port we're looking for.
                    if(NPort.ContainsPort(el, port))
                        boolean = 1;
                        return;
                    end
                elseif(isa(el, 'TerminatedTLine'))
                    % TerminatedTLine should not exist due to being flattened into a TLine.
                    error('TerminatedTLine should be flattened into TLine.');
%                     % Check if the transmission line or the terminator contains the port we're looking for.
%                     if(el.terminator == port || NPort.ContainsPort(el.tline, port))
%                         boolean = 1;
%                         return;
%                     end
                end
            end
        end
    end
end







































