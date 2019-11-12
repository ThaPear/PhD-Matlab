% classdef Line_Trenched < Line_Trenched_Dualpol_Bowtie
classdef Line_Trenched < Line_Trenched_Bowtie
    methods
        function this = Line_Trenched(er, L, erreal)
%             this@Line_Trenched_Dualpol_Bowtie(er, L, erreal);
            this@Line_Trenched_Bowtie(er, L, erreal);
        end
    end
end