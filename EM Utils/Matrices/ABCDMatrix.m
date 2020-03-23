classdef ABCDMatrix
    properties
        A
        B
        C
        D
    end
    methods
        function this = ABCDMatrix(A, B, C, D)
            switch(nargin)
                case 1 % The first argument is an ABCD matrix.
                    this.A = A.A;
                    this.B = A.B;
                    this.C = A.C;
                    this.D = A.D;
                case 4 % Simply 4 specified paremeters.
                    this.A = A;
                    this.B = B;
                    this.C = C;
                    this.D = D;
                otherwise
                    error('%s::ABCDMatrix:\n\tInvalid number of constructor arguments.', mfilename);
            end
        end
        
        function [result] = mul(this, ABCD)
            % Places the given matrix in series with this object.
            
            %        |------|    |------|
            %     -- |      | -- |      | -- 
            % -->    | this |    | ABCD |    -->
            %     -- |      | -- |      | -- 
            %        |------|    |------|
            
            % |--------------|
            % | AA+BC  AB+BD |
            % | CA+DC  CB+DD |
            % |--------------|
            
%             Anew = this.A .* ABCD.A + this.B .* ABCD.C;
%             Bnew = this.A .* ABCD.B + this.B .* ABCD.D;
%             Cnew = this.C .* ABCD.A + this.D .* ABCD.C;
%             Dnew = this.C .* ABCD.B + this.D .* ABCD.D;
%             result = ABCDMatrix(Anew, Bnew, Cnew, Dnew);
            
            result = ABCDMatrix(this.A .* ABCD.A + this.B .* ABCD.C, ...
                                this.A .* ABCD.B + this.B .* ABCD.D, ...
                                this.C .* ABCD.A + this.D .* ABCD.C, ...
                                this.C .* ABCD.B + this.D .* ABCD.D);
        end
        
%         function [result] = parallel(this, ABCD)
%             % Places the given matrix in parallel with this object.
%             
%             %              |------|
%             %          /-- |      | --\
%             %         /    | this |    \
%             %        / /-- |      | --\ \
%             %     --/ /    |------|    \ \--
%             % -->    X                  X    -->
%             %     --\ \    |------|    / /--
%             %        \ \-- |      | --/ /
%             %         \    | ABCD |    /
%             %          \-- |      | --/
%             %              |------|
%             
%             
%         end
    end
end