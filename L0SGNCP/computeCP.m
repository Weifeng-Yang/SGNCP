%% Calculate the value of the objective function
function [loss,S]=computeCP(var,ngmar)
X=ktensor(var);
S=tensor(X);
loss=norm(ngmar-S);
end

