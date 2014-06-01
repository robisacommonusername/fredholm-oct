%Return function S which can be evaluated anywhere. Specify a function chi
%and kernel H. Assume domain of integration is [0,1]
function S = apply_operator(H,chi)
	S = @(k) quad(@(z) H(k,z) * chi(z), 0, 1);
end
