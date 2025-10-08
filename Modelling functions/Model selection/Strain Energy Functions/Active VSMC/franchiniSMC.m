function [sigmatt_rr,sigmazz_rr,W,Ctttt,Czzzz] = franchiniSMC(xV,lambdat,lambdaz,lambdar)
    
    E1 = 2*(xV(2)*(xV(5)*lambdat+xV(6)*lambdaz-2)+(1-2*xV(2))*(xV(5)*lambdat-1));
    sigmatt_rr = 2*xV(1)*(1-xV(2))*xV(5)*lambdat.*(1+xV(3)*E1-xV(4)*E1.^2);
    sigmazz_rr = 2*xV(1)*xV(2)*xV(6)*lambdaz.*(1+xV(3)*E1-xV(4)*E1.^2);

    W = 0;
    Ctttt = 0;
    Czzzz = 0;
end  