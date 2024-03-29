% Terminated transmission line. Is terminated with the given terminator.
%   terminator is expected to be an element with a valid GetInputImpedance.
%
%      ---[Z]--- <-- Z = terminator
%      |       |
%    -------------
%    |           |
%    |   tline   |
%    |           |
%    -------------
%      |       |
%      o       o
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

            % Handle identity ABCD matrices.
            ind = ABCD.isidentity();
            Zin(ind) = ZL(ind);
            
%             if(any(isnan(Zin)))
%                 breakpoint;
%             end
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
        
        function BuildCST(this, project, parentcomponent)
            this.tline.BuildCST(project, parentcomponent);
            this.terminator.BuildCST(project, parentcomponent);
        end
        function newline = Flatten(this)
            % Flattens into a TLine with the terminator in Shunt at the end.
            newline = TLine({this.tline, Shunt(this.terminator)});
            newline = newline.Flatten();
        end
        function flippedline = Flip(this)
            % A flipped TerminatedTLine will have the termination at the input, so place it there as
            % a shunt.
            flippedline = this.tline.Flip();
            flippedline.elements = [{Shunt(this.terminator)}, flippedline.elements];
        end
    end
end