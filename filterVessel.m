%% Vessel Filter: function description
function vesselness = filterVessel(volume)
	volume = (1 - volume)*255;
	if ~isa(volume,'double')
		volume = single(volume);
	end

	sigma = 0;
	alpha = 0.5;
	beta = 0.5;
	FrangiC = 500;

	% Calculate 3D hessian
	[Dxx,Dyy,Dzz,Dxy,Dxz,Dyz] = Hessian3D(volume,sigma);

	% Calculate eigen values
	[lambda1,lambda2,lambda3] = eig3volume(Dxx,Dyy,Dzz,Dxy,Dxz,Dyz);

	clear Dxx Dyy Dzz Dxy Dxz Dyz;

	% Calculate absolute values of eigen values
	lambdaAbs1=abs(lambda1);
	lambdaAbs2=abs(lambda2);
	lambdaAbs3=abs(lambda3);

	% The Vesselness Features
	Ra = lambdaAbs2./lambdaAbs3;
	Rb = lambdaAbs1./sqrt(lambdaAbs2.*lambdaAbs3);

	% Second order structureness. 
	S = sqrt(lambdaAbs1.^2 + lambdaAbs2.^2 + lambdaAbs3.^2);
	A = 2*alpha^2;
	B = 2*beta^2;
	C = 2*FrangiC^2;

	expRa = (1-exp(-(Ra.^2./A)));
	expRb =    exp(-(Rb.^2./B));
	expS  = (1-exp(-(S.^2./C)));

	clear S A B C Ra Rb;

	vesselness = expRa.*expRb.*expS;

	clear expRa expRb expS;

	vesselness = vesselness.*and(lambda2 < 0,lambda3 < 0);
	vesselness(~isfinite(vesselness)) = 0;
end

%% Hessian3D: function description
function [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(volume,sigma)
	if sigma > 0
		volume = imgaussian(volume,sigma);
	else
		volume = volume;
	end

	% Create first and second order diferentiations
	[Dx,Dy,Dz] = gradient(volume);

	[Dxx,Dxy,Dxz] = gradient(Dx);
	clear Dx;

	[~,Dyy,Dyz] = gradient(Dy);
	clear Dy;

	[~,~,Dzz] = gradient(Dz);
	clear Dz;
end