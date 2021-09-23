% A port with a given impedance. 
%      o       o
%      |       |
%      |       |
%      --o Z o-- < Z = impedance

classdef Port
    properties
        number
        impedance
    end
    methods(Static)
        function port = New(impedance)
            % Gives a new port with the given impedance without having to worry about numbering.
            persistent number;
            if(isempty(number))
                number = 1;
            end
            port = Port(number, impedance);
            number = number + 1;
        end
    end
    methods
        function this = Port(number, impedance)
            this.number = number;
            this.impedance = impedance;
        end
        function zin = GetInputImpedance(this, isTE, f, k0, kr)
            zin = this.impedance;
        end
% Should not be necessary since ports are always in Shunt or as termination of a line.
%         function ABCD = GetABCD(this, isTE, f, k0, kr)
%             shunt = Shunt(Impedance(this.impedance));
%             ABCD = shunt.GetABCD(isTE, f, k0, kr);
%         end
        function boolean = eq(port1, port2)
            boolean = (port1.number == port2.number);
            if(boolean && (port1.impedance ~= port2.impedance))
                error('Ports have equal number but different impedance.\n');
            end
        end
        function flippedport = Flip(this)
            flippedport = this;
        end
    end
end