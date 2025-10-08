function [sigmatt_rr,sigmazz_rr, W, Ctttt, Czzzz]=Zulliger2fibre(parameterV,lambdatV,lambdazV,lambdarV,activator)
% The function returns the circumferential and axial stresses, the stored
% elastic energy, and the small-on-large circumferential and axial material
% stiffnesses (NOT implemented yet) for Zulliger's constitutive model given a parameter
% vector parameterV, the deformation in the three principal directions, an activator 
% that activates/deactivates the contribution from each wal constituent.
%
% For a detailed desription of the constitutive model refer to Zulliger and Stergiopulos 
% Structural strain energy function applied to te aging of the human aorta  
% Journal of Biomechanics, 40, pp.3061-3069, 2007. Note that an additional
% collagen fibre family has been added.
%
% parameterV is the 12-element vector of model parameters where:
%   position 1: elastin stiffness-like parameter
%   position 2: collagen stiffness-like parameter
%   position 3: collagen stiffness-like parameter for compression
%   position 4: collagen fibre engagment distribution scale parameter (b)
%   position 5: collagen fibre engagment distribution shape parameter (k)
%   position 6: collagen fibre orientation (in radians)
%   position 7: elastin circumferential deposition stretch 
%   position 8: elastin axial deposition stretch 
%   position 9: collagen stiffness-like parameter
%   position 10: collagen stiffness-like parameter for compression
%   position 11: collagen fibre engagment distribution scale parameter (b)
%   position 12: collagen fibre engagment distribution shape parameter (k)
%
% For a detailed description of the small-on-large formulation refer to
% Baek et al. Theory of small on large: Potential utility in computations
% of fluid-solid interactions in arteries. Computer methods in applied
% mechanics nad engineering, 196(31-32), pp.3070-3078.

if(~exist('activator'))
    activator=[1,1]; % selective activation of collagen and elastin --> if not given as input, is set to [1, 1]
end

E_fibre_dir = 0.5*((lambdatV.^2*sin(parameterV(6))^2+lambdazV.^2*cos(parameterV(6))^2)-1); %strain in the fbre direction
I1 = (parameterV(7)*lambdatV).^2+(parameterV(8)*lambdazV).^2+(lambdarV/parameterV(7)/parameterV(8)).^2;

% Building the matrix of arguments of the probability distribution
% function for diagonal fibres
y_vect=E_fibre_dir+1;
x_vect=0:0.001:max(y_vect);
matB=y_vect-x_vect;
locator=find(matB<0);
matB(locator)=0;

% Building the matrix of arguments of the probability distribution
% function for circumferential fibres
y_vect=0.5*(lambdatV.^2-1)+1;
x_vect=0:0.001:max(y_vect);
matC=y_vect-x_vect;
locator=find(matC<0);
matC(locator)=0;

% Probability distribution function matrix for convolution
p_diagonal_fibre=parameterV(5)*parameterV(4)^parameterV(5)*(matB).^(parameterV(5)-1)./(parameterV(4)^parameterV(5)+(matB).^parameterV(5)).^2;
p_circ_fibre=parameterV(12)*parameterV(11)^parameterV(12)*(matC).^(parameterV(12)-1)./(parameterV(11)^parameterV(12)+(matC).^parameterV(12)).^2;

% Building the matrix of argument of the diagonal fibre stress function
q_vect=-1:0.001:max(E_fibre_dir);
matA=zeros(length(E_fibre_dir),1)+q_vect;
sigmatt_diagonal_fibre=(parameterV(2)*(matA>0)+parameterV(3)*(matA<=0)).*(2*matA)*sin(parameterV(6))^2;
sigmazz_diagonal_fibre=(parameterV(2)*(matA>0)+parameterV(3)*(matA<=0)).*(2*matA)*cos(parameterV(6))^2;

% Building the matrix of argument of the diagonal fibre stress function
s_vect=-1:0.001:max(0.5*(lambdatV.^2-1));
matD=zeros(length(lambdatV),1)+s_vect;
sigmatt_circ_fibre=(parameterV(9)*(matD>0)+parameterV(10)*(matD<=0)).*(2*matD);

% Calculate stress as sum of collagen and elastin
sigmatt_rr=activator(1)*parameterV(1)*((lambdatV*parameterV(7)).^2-(lambdarV/parameterV(7)/parameterV(8)).^2)+activator(2)*2*lambdatV.^2.*...
    (sum((sigmatt_diagonal_fibre.*p_diagonal_fibre*0.001)')'+sum((sigmatt_circ_fibre.*p_circ_fibre*0.001)')');
sigmazz_rr=activator(1)*parameterV(1)*((lambdazV*parameterV(8)).^2-(lambdarV/parameterV(7)/parameterV(8)).^2)+activator(2)*2*lambdazV.^2.*...
    sum((sigmazz_diagonal_fibre.*p_diagonal_fibre*0.001)')';

% Calculate stored energy
W_diagonal_fibre=0.5*(parameterV(2)*(matA>0)+parameterV(3)*(matA<=0)).*matA.^2;
W_circ_fibre=0.5*(parameterV(2)*(matD>0)+parameterV(3)*(matD<=0)).*matD.^2;
W=activator(1)*(parameterV(1)/2)*(I1-3)+activator(2)*(sum((W_diagonal_fibre.*p_diagonal_fibre*0.001)')+sum((W_circ_fibre.*p_circ_fibre*0.001)'));

% Calculate the stiffness
Ctttt=2*(sigmatt_rr+activator(1)*parameterV(1)*((lambdarV/parameterV(7)/parameterV(8)).^2)) +...
    activator(2)*4*sin(parameterV(6))^4*lambdatV.^4.*sum(((parameterV(2)*(matA>0)+parameterV(3)*(matA<=0)).*p_diagonal_fibre*0.001)');
Czzzz=2*(sigmazz_rr+activator(1)*parameterV(1)*((lambdarV/parameterV(7)/parameterV(8)).^2)) +...
    activator(2)*4*cos(parameterV(6))^4*lambdazV.^4.*sum(((parameterV(2)*(matA>0)+parameterV(3)*(matA<=0)).*p_diagonal_fibre*0.001)');

end