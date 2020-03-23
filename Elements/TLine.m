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
                
                element = this.elements{i};
                ABCDelement = element.GetABCD(isTE, f, k0, kr);
                if(max(isnan(ABCDelement.A(:))) || max(isnan(ABCDelement.B(:))) || max(isnan(ABCDelement.C(:))) || max(isnan(ABCDelement.D(:))))
                    breakpoint;
                    element.GetABCD(isTE, f, k0, kr);
                end
                % Multiply the ABCD for each element together.
                ABCD = ABCD.mul(ABCDelement);
                if(max(isnan(ABCD.A(:))) || max(isnan(ABCD.B(:))) || max(isnan(ABCD.C(:))) || max(isnan(ABCD.D(:))))
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
%             
%             etaTE  = sqrt(((1+Ste.s11 ).^2 - Ste.s21.^2 ) ./ ((1-Ste.s11 ).^2 - Ste.s21.^2 )) .* sec(th );
%             etaTE0 = sqrt(((1+Ste0.s11).^2 - Ste0.s21.^2) ./ ((1-Ste0.s11).^2 - Ste0.s21.^2)) .* sec(th0);
%             
%             etaTM  = sqrt(((1+Stm.s11 ).^2 - Stm.s21.^2 ) ./ ((1-Stm.s11 ).^2 - Stm.s21.^2 )) .* cos(th );
%             etaTM0 = sqrt(((1+Stm0.s11).^2 - Stm0.s21.^2) ./ ((1-Stm0.s11).^2 - Stm0.s21.^2)) .* cos(th0);
%             
%             if(real(etaTE)  < 0); etaTE  = -etaTE;  end
%             if(real(etaTM)  < 0); etaTM  = -etaTM;  end
%             if(real(etaTE0) < 0); etaTE0 = -etaTE0; end
%             if(real(etaTM0) < 0); etaTM0 = -etaTM0; end
%             
%             zetaTE  = Ste.s21  ./ (1 - Ste.s11 .*((etaTE .*cos(th ) - 1) ./ (etaTE  .* cos(th ) + 1)));
%             zetaTE0 = Ste0.s21 ./ (1 - Ste0.s11.*((etaTE0.*cos(th0) - 1) ./ (etaTE0 .* cos(th0) + 1)));
%             
%             zetaTM  = Stm.s21  ./ (1 - Stm.s11 .*((etaTM ./cos(th ) - 1) ./ (etaTM  ./ cos(th ) + 1)));
%             zetaTM0 = Stm0.s21 ./ (1 - Stm0.s11.*((etaTM0./cos(th0) - 1) ./ (etaTM0 ./ cos(th0) + 1)));
%             
%             nTE  = @(m) (sqrt( ( (log(abs(zetaTE )) + 1j.*(angle(zetaTE ) + 2*pi*m)) ./ (-1j.*k0.*d)).^2 + sin(th ).^2));
%             nTE0 = @(m) sqrt( ( (log(abs(zetaTE0)) + 1j.*(angle(zetaTE0) + 2*pi*m)) ./ (-1j.*k00.*d)).^2 + sin(th0).^2);
%             
%             nTM  = @(m) sqrt(((log(abs(zetaTM )) + 1j.*(angle(zetaTM ) + 2*pi*m)) ./ (-1j.*k0.*d)).^2 + sin(th ).^2);
%             nTM0 = @(m) sqrt(((log(abs(zetaTM0)) + 1j.*(angle(zetaTM0) + 2*pi*m)) ./ (-1j.*k00.*d)).^2 + sin(th0).^2);
%             
%             
% %             eta  = @(S, th)     sqrt(((1+S.s11).^2 - S.s21.^2) ./ ((1-S.s11).^2 - S.s21.^2)) .* sec(th);
% %             zeta = @(S, th)     S.s21 ./ (1 - S.s11.*((eta(S, th).*cos(th) - 1) ./ (eta(S, th) .* cos(th) + 1)));
% %             n    = @(S, th, m)  sqrt(((log(abs(zeta(S, th))) + 1j.*angle(zeta(S, th)) + 2*pi*m) ./ (-1j.*k0.*d)).^2) + sin(th).^2;
%             
%             % Epsilon as function of m.
%             erx  = @(m)         nTM0(m) ./ etaTM0;
%             ery  = @(m)         nTE0(m) ./ etaTE0;
%             erz  = @(m)         erx(m) .* sin(th).^2 ./ (sin(th).^2 - nTM(m).^2 + nTM0(m).^2);
%             erz2 = @(m)         sin(th).^2 ./ (nTM0(m) .* etaTM0 - erx(m) .* etaTM.^2);
%             
%             mux  = @(m)         (nTE0(m) .* etaTE0);
%             muy  = @(m)         nTM0(m) .* etaTM0;
%             % ERROR
%             muz  = @(m)         sin(th).^2 ./ ((nTE0(m) ./ etaTE0) - mux(m) ./ etaTE.^2);
%             muz2 = @(m)         mux(m) .* sin(th).^2 ./ (sin(th).^2 - nTE(m).^2 + nTE0(m).^2);
%             
%             printvar muz(-3:3);
%             printvar muz2(-3:3);
%             ms = -3:1:3;
%             
%             errors = abs(erz(ms) - erz2(ms));
%             [~, idx] = min(errors); % Get index of lowest value
%             mbar = ms(idx);         % Get corresponding m.
%             
%             errors = abs(muz(ms) - muz2(ms));
%             [~, idx] = min(errors); % Get index of lowest value
%             m = ms(idx);            % Get corresponding m.
%             
% %             printvar etaTE;
% %             printvar etaTE0;
% %             figureex(3);
% %                 plot(real(nTE0(ms)));
% %                 plot(imag(nTE0(ms)), '--');
% %                 plot(real(mux(ms)));
% %                 plot(imag(mux(ms)), '--');
% %                 legend('nTE0', 'nTE0im', 'mux', 'muxim');
% %             printvar muz(-3:3)
% %             printvar muz2(-3:3)
%             
% 
%             % The error minimization doesn't work correctly.
%             % Use m = mbar = 0 for all calculations.
% %             mbar = 0;
% %             m = 0;
% 
%             er.x = erx(mbar);
%             er.y = ery(mbar);
%             er.z = erz(mbar);
%             mu.x = mux(m);
%             mu.y = muy(m);
%             mu.z = muz2(m);
            
        end
        function PlotEpsilons(this, fs, hAx)
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
        function BuildCST(this, project)
            for(i = 1:length(this.elements))
                this.elements{i}.BuildCST(project);
            end
        end
        
        % Allow indexing of elements by directly doing tline{3} for
        % tline.elements{3}
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