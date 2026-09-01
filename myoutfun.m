function [stop,options,optchanged] = myoutfun(x,optimValues,state)

persistent Xhist Fhist

stop = false;
optchanged = false;

switch state
    case 'init'
        Xhist = [];
        Fhist = [];

    case 'iter'
        Xhist = [Xhist; x(:)'];
        Fhist = [Fhist; optimValues.fval];

    case 'done'
        assignin('base','Xhist',Xhist);
        assignin('base','Fhist',Fhist);
end
end