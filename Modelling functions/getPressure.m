function P = getPressure(parameterV_elastic,ri,lambdazV,h_ref,ro_ref,l_ref,fun)  
    ri_ref = ro_ref-h_ref;
    vessel_volume = pi*(ro_ref^2-ri_ref^2)*l_ref;

    ro = radiusFromIincompressibility(ri,lambdazV*l_ref,vessel_volume,2);
    
    h = ro-ri;
    lambdatV = (ri+h/2)./(ri_ref+h_ref/2);
    lambdarV = 1./(lambdazV.*lambdatV);

    [sigmatt_rrV,sigmazz_rrV,~,~,~] = fun(parameterV_elastic,lambdatV,lambdazV,lambdarV);
    [P,~] = stress2pressure_force(sigmatt_rrV,sigmazz_rrV,(ri+h/2),h);
end