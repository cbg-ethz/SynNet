% this is basically the matlbas bionomial random generator, except that it
% uses normal and binomial approximation for extreme cases.
% NOTE this is written only for the case that there's a matrix N and a
% matrix P and the result is aa matrix R with the same size.

% Pejman May 2015
function R = Pej_Binornd(N, P)
F1 = N<1000;
F2 = P>.05 & P<.95;

R = nan(size(N));

R(F1) = binornd(N(F1), P(F1));

F3 = ~F1 & F2;
npF3 = N(F3).*P(F3);
R(F3) = round(randn(sum(F3(:)), 1).* sqrt(npF3.*(1-P(F3))) + npF3); % normal approximation

F4 = ~F1 & ~F2 & P<=.05;
npF4 = N(F4).*P(F4);
R(F4)= poissrnd(npF4); % poisson approximation

F5 = ~F1 & ~F2 & P>=.95;
npF5 = N(F5).*(1-P(F5));
R(F5)= N(F5) - poissrnd(npF5); % poisson approximation

Ffail = R<0 | R>N; % these cases can happen in approximates, I exclude them and re-do them, this makes it a truncated distribution
if any(Ffail(:))
    R(Ffail) = Pej_Binornd(N(Ffail), P(Ffail));
end
end