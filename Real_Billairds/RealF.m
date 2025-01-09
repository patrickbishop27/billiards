% Iterative Billiard map for real value case. Takes in current theta value
% and r value and outputs next r value and theta value.  
% INPUT: input = [r, theta] where r = cos(gamma) and theta is parameter value 
%                corresponging to initial point on boundary.
%        coefx = the coefficients for the cosine terms in the x coordinate 
%        for the table boundary.
%        coefy = coefficients for the sine terms in the y coordinate for
%        the table boundary.
% OUTPUT: [rhat, thetahat] for the next point of contact of the billiard ball 
%         with the boundary.
function output = RealF(input,coefx,coefy)

th1 = input(2);
r1   = input(1);
maxresid = 0;
maxk = 0;
flag = 0; 

[~,dz1] = Bill_Table(th1,coefx,coefy);
v1 = r2v2(r1,th1,dz1/norm(dz1));
v=priorv(v1,dz1);

% Create guess for next newtons method
[distguess,thguess] = nextstepell(th1,v,coefx(1),coefy(1));	
guess = [distguess;thguess];    
% Construct function from vector v and boundary. Use Newton's method to
% find point of intersection.
[pts2,v,~] = findbouncevector(th1,v,coefx,coefy);
L =@(s) pts2 + s.*v;

F =@(x) L(x(1)) - Bill_Table(x(2),coefx,coefy);
dF = @(x) [v,-dBill_Table(x(2),coefx,coefy)];

% Apply deflation to avoid converging to point of origin.     
deflat =  @(x) (1/norm([x(1)-0;x(2)-th1])^2);
ddeflat = @(x) -2*[x(1),x(2)-th1]*deflat(x)^2;
	
G = @(x) F(x)*deflat(x);
dG =@(x) F(x)*ddeflat(x)+deflat(x)*dF(x);

[newtval,~,kk] = Newtons(G,dG,guess);

% Find new r value and save both [r,th]
th2 = newtval(2);
[~,dz2] = Bill_Table(th2,coefx,coefy);	
r2 = v2r(v,dz2); 
output = [r2;th2];


%% nextstepell
% Finds the next point of contact on the ellipse. 
function [tell,th2] = nextstepell(th1,v1,a,b)
	%Use the previous v so we can figure out the ellipse bounce instead of the actual bounce
	tpi = 2*pi;
	
	[~,vell] = findbouncevector(th1,v1,a,b);

	dell = vell./[a;b];
	tell = -2 * (dell(1)*cos(tpi*th1)+dell(2)*sin(tpi*th1))/(dell(1)^2+dell(2)^2);
	th2 = mod(atan2(sin(tpi*th1)+tell*dell(2),cos(tpi*th1)+tell*dell(1)),tpi)/tpi;

	
%% findbouncevector
% Law of reflection
function [pts,v,r] = findbouncevector(th,v1,coefx,coefy)
    % Find normal vector to p2
    [pts,n1] = Bill_Table(th,coefx,coefy);
    n1 = n1/norm(n1);
    n = [-n1(2);n1(1)];
    v1 = v1/norm(v1);
    r = dot(v1,n1);
    v2 = v1-2*dot(n,v1)*n;
    v = v2/norm(v2);
