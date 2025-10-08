function [sigmatt_rr,sigmazz_rr,W,Ctttt,Czzzz] = ZulligerSMCdispersion(xV,lambdat,lambdaz,lambdar)

    sigmatt_rr = xV(1)*(1-xV(2))*(lambdat*xV(3)-1)+xV(1)*xV(2)*(xV(3)*lambdat)./(xV(3)*lambdat+xV(4)*lambdaz);
    sigmazz_rr = xV(1)*xV(2)*lambdaz*xV(4).*(1-1./(xV(3)*lambdat+xV(4)*lambdaz));

    W = xV(1)*(1-2*xV(2))*(xV(3)*lambdat+log(xV(3)*lambdat)-1)+xV(1)*xV(2)*((xV(3)*lambdat+xV(4)*lambdaz)-...
        log((xV(3)*lambdat+xV(4)*lambdaz)/2)-2);

    Ctttt = xV(1)*(1-2*xV(2)+xV(2)*(xV(3)*lambdat).^2./(xV(3)*lambdat+xV(4)*lambdaz).^2);
    Czzzz = xV(1)*xV(2)*(xV(4)*lambdaz).^2./(xV(3)*lambdat+xV(4)*lambdaz).^2;
end