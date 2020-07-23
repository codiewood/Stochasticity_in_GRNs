% syms ha hr fr fa positive;
inputs = {1e-4,1e-2,2,1e-1};

%sets defaults if not all args in input
numargs = length(inputs);
args = {1e-4,1e-2,2,1e-1};
args(1:numargs) = inputs;
[hr,fr,ha,fa] = args{:};

M = [-(ha+hr) ha 0 hr; fa -(fa+hr) hr 0;0 fr -(fr+fa) fa;fr 0 ha -(ha+fr)]
[V,D] = eig(M');
E = diag(D)
[A,I] = min(abs(E));
s = V(:,I);
n = sum(s);
st = s/n

