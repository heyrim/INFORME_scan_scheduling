function f = logfun(x)

global FHIST XHIST

f = originalObjective(x);

FHIST(end+1,1) = f;
XHIST(end+1,:) = x(:)';
end