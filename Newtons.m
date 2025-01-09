%Multivariable Newton's Method
function [xfin,stop,k] = Newtons(F,dF,guess)   

    x = guess;
    stop = 1;
    k = 1;

    while (stop > 10^(-14)) && (k<1000)
		diff = -dF(x)\F(x);
        x = x + diff; 
        stop = norm(diff,'inf');
        k = k+1;
     end
    
	extraits = 2;
    for j = 1:extraits
        diff = -1*dF(x)\F(x);
        x = x + diff; 
        stop = norm(diff,'inf');
        k = k+1;
     end


        
     xfin = x;
     
end