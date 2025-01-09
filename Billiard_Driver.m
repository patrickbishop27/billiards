%% Driver for analytic billiard map. 
% Input: numvals = number of r values taken at initial point on the
%        boundary. r = cos(gamma) where gamma is the angle made between the 
%        the initial trajectory of the ball and the tangent line at point
%        of origin. 
%        numits = number of iterations in each orbit.
%        thetainit = initial point on the boundary where theta is between 0
%        and 1. 
%        coefx = the coefficients for the cosine terms in the x coordinate 
%        for the table.
%        coefy = coefficients for the sine terms in the y coordinate for
%        the table.
% Output: Phase portrait for the billiard table colored by rotation number
% and dig value. 
% Example: Billiard_Driver(50,100,0.23,[1.1,0.03],[1,0.028]) will plot the phase space for 50 orbits 
% each with length 100 all starting from theta = 0.23 on the table with
% equation x = 1.1*cos(2*pi*t) + 0.03*cos(2*pi*t)
%          y = 1*sin(2*pi*t) + 0.028*sin(2*pi*t) 

function Billiard_Driver(numvals,numits,thetainit,coefx,coefy)


rvec = linspace(-0.99,0.99,numvals); 
rfinal = zeros(numvals,numits+1);
thfinal= zeros(numvals,numits+1);
rotfinal = zeros(numvals,numits+1);
digfinal = zeros(numvals,numits+1);


for j = 1:numvals
    rit = zeros(numits,1);
    thetait = zeros(numits,1);
    rit(1) = rvec(j);
	thetait(1) = thetainit; 

	for k = 1:numits
		output = RealF([rit(k),thetait(k)],coefx,coefy);
		rit(k+1) = output(1);
		thetait(k+1) = mod(output(2),1);
    end

    thfinal(j,:) = thetait(1:end)';
    rfinal(j,:) = rit(1:end)';

    [rotn,diggn] = Extra_info2(thfinal(j,1:end-1),numits);

    digfinal(j,:) = ((1-isnan(diggn))*diggn)*ones(1,numits+1);
	rotfinal(j,:) = min(mod(rotn,1),1)*ones(1,numits+1);


end


rplot = reshape(rfinal,(numits+1)*numvals,1);
thplot = reshape(thfinal,(numits+1)*numvals,1);
rotplot = reshape(rotfinal,(numits+1)*numvals,1);
digplot = reshape(digfinal,(numits+1)*numvals,1);
thplot2 = mod(thplot,1);

%% Rotation Number
%figure; 
%scatter(thplot2,rplot,1,rotplot);
%caxis([0,1]);
%colormap(hsv);
%colorbar;
%xlabel('\fontsize{30} \theta')
%ylabel('\fontsize{30} r')


%% Dig Value
figure; 
scatter(thplot2,rplot,1,digplot);
colormap(flipud(jet));
caxis([0,12])
colorbar;
%xlabel('\fontsize{25} \theta')
%ylabel('\fontsize{25} r')

%% Rotation without chaos
trust = find(digplot>4.875);
dontrust = find(digplot<=4.875);
figure; 
scatter(thplot2(trust),rplot(trust),1,rotplot(trust)); hold on
scatter(thplot2(dontrust),rplot(dontrust),1,'k');
caxis([0,1])
colormap(hsv)
colorbar;
%xlabel('\fontsize{20} \theta')
%ylabel('\fontsize{20} r')
%% Histogram of Dig
%figure;
%histogram(digplot,30,"Normalization","probability");
%xlabel('\fontsize{20} dig_{T}')
%ylabel('\fontsize{20} Proportion')
