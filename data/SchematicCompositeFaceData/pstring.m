function [pstr, ptag] = pstring(p)

if p <= .1 && p > .05
    tp   = sprintf('%2.2f',  p);
    [a, b] = strtok(tp, '.');
    pstr = ['= ' b];
    ptag = '+-';
elseif p <= .05 && p > .01
    pstr = '< .05';
    ptag = '*';
elseif p <= .01 && p > .001
    pstr = '< .01';
    ptag = '**';
elseif p <= .001 
    pstr = '< .001';
    ptag = '***';
else
    pstr = '> .05';
    ptag = '';
end