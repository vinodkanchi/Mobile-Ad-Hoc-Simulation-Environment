function m = f_msequence_gen(start, coeff, len)
% f_msequence_gen(starting value, coefficient ,length)
%
% The function f_msequence_gen(starting value, coefficient ,length)
% generates cyclic codes and m-sequences. 
%
% input variables
% - starting value is the beginning sequence of the cyclic code x^n --> x^1
% - coefficient represents the coefficient of the generating polynomial
% - the length means the number of the generated bits
%
% output variable
% - the output is a vector that contains the cyclic code
%
% example: m1=f_msequence_gen([0 0 0 0 1], [1 0 0 1 0 1], 31)
  
if (length(start)+1) == length(coeff)
   
    % init
    start = fliplr(start);
    coeff = fliplr(coeff);
    m = zeros(1, len);

    for k= 0 : len-1
        % m is the 5. value of start
        m(k+1) = start(length(start));

        % binary
        sc=start*coeff(2:length(coeff))';
        % addition
        first=mod(sc,2);
        
        start=[first, start(1:length(start)-1)];
    end %k
else 
   error('The length of the starting value has to be one more than the length of the coefficient of the polynomial');
end   