% Terminated transmission line. Is terminated with the given terminator.
%   terminator is expected to be an element with a valid GetInputImpedance.

classdef TerminatedTLine < Element
    properties
        tline
        terminator
    end
    methods
        function this = TerminatedTLine(tline, terminator)
            this.tline = tline;
            this.terminator = terminator;
        end
        
        function Zin = GetInputImpedance(this, isTE, f, k0, kr)
            ABCD = this.tline.GetABCD(isTE, f, k0, kr);
            
            Zmat = ABCD2Z(ABCD);
            
            
            ZL = this.terminator.GetInputImpedance(isTE, f, k0, kr);
            
            % Convert Z-matrix with load into input impedance.
            Zin = Zmat.z11 - (Zmat.z12 .* Zmat.z21) ./ (Zmat.z22 + ZL);
%             zout = Zmat.z22 - (Zmat.z12 .* Zmat.z21) ./ (Zmat.z11 + ZL);
            if(max(isnan(Zin)))
                breakpoint;
            end
        end
        
        function [ABCD] = GetABCD(this, isTE, f, k0, kr)
            error('%s::GetABCD:\n\tShould not be called on TerminatedTLine.\n\tABCD matrix is not valid.', mfilename);
        end
        
        function h = GetHeight(this)
            h = 0;
            h = h + this.tline.GetHeight();
            h = h + this.terminator.GetHeight();
        end
        
        function h = GetEffectiveHeight(this, f)
            h = 0;
            h = h + this.tline.GetEffectiveHeight(f);
            h = h + this.terminator.GetEffectiveHeight(f);
        end
        
        function BuildCST(this, project)
            this.tline.BuildCST(project);
            this.terminator.BuildCST(project);
        end
    end
end