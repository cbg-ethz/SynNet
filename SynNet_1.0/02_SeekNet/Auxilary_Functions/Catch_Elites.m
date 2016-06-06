% This finds the Nth unique highest number in the performance
% array, and returns an array for those larger than it.

function [Elite_Idx, nElites] = Catch_Elites(Performance, N)
        [Sp, SpI]  = sort(Performance, 'descend');
        dSp = Sp(1:(end-1)) - Sp(2:end);
        Jumps = [find(dSp>0); length(SpI)];
        nElites = Jumps(min(N, length(Jumps)));
        Elite_Idx =  SpI(1:nElites);
end