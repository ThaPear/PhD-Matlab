classdef TLine < Element
    properties
        elements
    end
    properties(Dependent)
        N % Number of elements in line.
    end
    methods
        function this = TLine(elements)
            if(nargin >= 1)
                this.elements = elements;
            end
        end
        function n = get.N(this)
            n = length(this.elements);
        end
        function Zin = GetInputImpedance(this, isTE, f, k0, kr)
            ZL = this.elements{end}.GetInputImpedance(isTE, f, k0, kr);
            
            line = TLine(this.elements(1:end-1));
            if(line.N > 0)
                ABCD = line.GetABCD(isTE, f, k0, kr);
            else
                Zin = ZL;
                return;
            end
            
            Zmat = ABCD2Z(ABCD);
            
            % Convert Z-matrix with load into input impedance.
            Zin = Zmat.z11 - (Zmat.z12 .* Zmat.z21) ./ (Zmat.z22 + ZL);
        end
        function [ABCD] = GetABCD(this, isTE, f, k0, kr)
            ABCD = ABCDMatrix(1, 0, ...
                              0, 1);
            for(i = 1:length(this.elements))
                % Cascade the ABCD matrices of the elements in the transmission line.
                element = this.elements{i};
                ABCDelement = element.GetABCD(isTE, f, k0, kr);
                if(any(isnan(ABCDelement.A(:))) || any(isnan(ABCDelement.B(:))) || any(isnan(ABCDelement.C(:))) || any(isnan(ABCDelement.D(:))))
                    breakpoint;
                    element.GetABCD(isTE, f, k0, kr);
                end
                % Multiply the ABCD for each element together.
                ABCD = ABCD.mul(ABCDelement);
                if(any(isnan(ABCD.A(:))) || any(isnan(ABCD.B(:))) || any(isnan(ABCD.C(:))) || any(isnan(ABCD.D(:))))
                    breakpoint;
                end
%                 dispex('%s\n', class(element));
%                 ABCD
            end
        end
        function h = GetHeight(this)
            h = 0;
            for(i = 1:length(this.elements))
                h = h + this.elements{i}.GetHeight();
            end
        end
        function h = GetEffectiveHeight(this, f)
            h = 0;
            for(i = 1:length(this.elements))
                h = h + this.elements{i}.GetEffectiveHeight(f);
            end
        end
        function [er, mu] = GetEpsilonMu(this, f, th, ph)
            th0 = eps;
            if(th <= eps)
                th = 10*pi/180;
            end
            
            [k0, ~, ~, ~] = k(f, 1, th, ph);
            [k00, ~, ~, ~] = k(f, 1, th0, ph);
            kr  = k0.*sin(th );
            kr0 = k0.*sin(th0);
            
            kz = -1j*sqrt(-(k0^2-kr^2));
            
            % Calculate the values for an angle of incidence different than 0.
            ABCDte = this.GetABCD(1, f, k0, kr);
            ABCDtm = this.GetABCD(0, f, k0, kr);
            
            [kd, ~, ~, kzd] = k(f, 1, th, ph);
            [~, z0te, z0tm] = z(1, kd, kzd);
            Ste = ABCD2S(ABCDte, z0te, z0te);
            Stm = ABCD2S(ABCDtm, z0tm, z0tm);
            
            % Calculate the values for an angle of incidence of 0.
            ABCDte0 = this.GetABCD(1, f, k00, kr0);
            ABCDtm0 = this.GetABCD(0, f, k00, kr0);
            Ste0 = ABCD2S(ABCDte0);
            Stm0 = ABCD2S(ABCDtm0);
            
            % Height of the structure under test.
            d = this.GetHeight();
            [er, mu] = EpsilonMu(f, Ste, Stm, Ste0, Stm0, d, th0, th, ph);
        end
        function PlotEpsilons(this, fs, hAx)
            % Plot the permittivity of this line as a function of the given frequencies.
            if(nargin < 3 || isempty(hAx))
                hFig = figureex;
                hAx = hFig.CurrentAxes;
                xlabel('Frequency [GHz]');
                ylabel([char(949), 'r']);
                title([char(949), 'r vs frequency']);
            end
            for(i = 1:length(this.elements))
                this.elements{i}.PlotEpsilons(fs, hAx);
            end
        end
        function BuildCST(this, project, parentcomponent)
            if(nargin < 3 || isempty(parentcomponent))
                parentcomponent = '';
            end
            % Build the sub-elements in CST.
            for(i = 1:length(this.elements))
                this.elements{i}.BuildCST(project, parentcomponent);
            end
        end
        function flippedtline = Flip(this)
            % Flip the transmission line upside-down.
            flippedtline = TLine(this.elements(length(this.elements):-1:1));
            flippedelements = this.elements(length(this.elements):-1:1);
            for(i = 1:length(flippedelements))
                flippedelements{i} = flippedelements{i}.Flip();
            end
        end
        function newline = Flatten(this)
            newline = this;
            restart = 1;
            while(restart)
                % Assume we don't need to restart.
                restart = 0;
                for(iel = 1:length(newline.elements))
                    el = newline.elements{iel};
                    if(strcmp(class(el), 'TLine')) %#ok<STISA> 'isa' also accepts subclasses, which we don't want.
                        % Flatten the sub-TLine.
                        el = el.Flatten();
                        % Insert the newly flattened TLine into this one.
                        newline.elements = [newline.elements(1:iel-1), el.elements, newline.elements(iel+1:end)];
                        % Restart the iteration of the tline.
                        restart = 1;
                        break;
                    elseif(isa(el, 'Shunt') || isa(el, 'Series'))
                        if(strcmp(class(el.impedance), 'TLine')) %#ok<STISA> 'isa' also accepts subclasses, which we don't want.
                            newline.elements{iel}.impedance = el.impedance.Flatten();
                        elseif(isa(el.impedance, 'TerminatedTLine'))
                            newline.elements{iel}.impedance = el.impedance.Flatten();
                        end
                    elseif(isa(el, 'TerminatedTLine'))
                        % The TerminatedTLine flattens into a TLine.
                        el = el.Flatten();
                        % Insert the newly flattened TLine into this one.
                        newline.elements = [newline.elements(1:iel-1), el.elements, newline.elements(iel+1:end)];
                        % Restart the iteration of the tline.
                        restart = 1;
                        break;
                    end
                end
            end
        end
        % Allow indexing of elements by directly doing tline{3} instead of tline.elements{3}
        % Disabled because it made indexing very slow.
%         function varargout = subsref(this, S)
%             if(S(1).type(1) == '{')%(strcmp(S(1).type, '{}'))
%                 if(length(S.subs) == 1)
%                     varargout{1:nargout} = this.elements{S.subs{1}};
%                 else
%                     el = this.elements{S.subs{1}};
%                     if(~isprop(el, 'elements'))
%                         error('Cannot index non-TLine element %s as cell', class(el));
%                     end
%                     varargout{1:nargout} = el.subsref(struct('type', '{}', 'subs', {S.subs(2:end)}));
%                 end
%                 return;
%             end
%             [varargout{1:nargout}] = builtin('subsref', this, S);
%         end
    end
end