function test
A = rand(40,40) + 1j*eye(40);
x = rand(40,1);
b = A*x;
x0 = zeros(40,1);
[x,resvec]=idrs(A,b,4,1e-6,100,x0);

resvec

end

function [x,resvec]=idrs(A,b,s,tol,maxit,x0);
%--------------- Creating start residual: ----------
N = length(b);
x = x0;
r = b - A*x;
normr = norm(r);
tolr = tol * norm(b);
resvec=[normr];% tol: relative tolerance
if (normr <= tolr)
  iter=0;
  return;
end;% Initial guess is a good enough solution
%----------------- Shadow space: --------------------
rand("state", 0);%for reproducibility reasons.
P = rand(N,s);
P(:,1) = r;
P = orth(P)';
% Only for comparison with Bi-CGSTAB
% transpose for efficiency reasons.
%---------------- Produce
dR = zeros(N,s);
dX = zeros(N,s);
for k = 1:s
  v = A*r;
  om = dot(v,r)/dot(v,v);
  dX(:,k) = om*r;
  dR(:,k) = -om*v;
  x = x + dX(:,k);
  r = r + dR(:,k);
  normr = norm(r);
  resvec = [resvec;normr];
  M(:,k) = P*dR(:,k);
end

iter = s;
oldest = 1;
m = P*r;
while ( normr > tolr ) & (iter < maxit )
  for k = 0:s
    c = M\m;
    q = -dR*c;			% s-1 updates + 1 scaling
    v = r + q;			% simple addition
    if ( k == 0 ) % 1 time:
      t = A*v;			% 1 matmul
      om = dot(t,v)/dot(t,t);	% 2 inner products
      dR(:,oldest) = q - om*t;  % 1 update
      dX(:,oldest) = -dX*c + om*v;% s updates + 1 scaling
    else
      dX(:,oldest) = -dX*c + om*v;% s updates + 1 scaling
      dR(:,oldest) = -A*dX(:,oldest); % 1 matmul
    end
    r = r + dR(:,oldest);
    x = x + dX(:,oldest);
    iter = iter + 1;% simple addition
    normr=norm(r);
    resvec = [resvec;normr];% 1 inner product (not counted)
    dm = P*dR(:,oldest);
    M(:,oldest) = dm;
    m = m + dm;% s inner products
	       % cycling s+1 times through matrices with s columns:
    oldest = oldest + 1;
    if ( oldest > s )
      oldest = 1;
    end
  end; % k = 0:s
end; %while

end
