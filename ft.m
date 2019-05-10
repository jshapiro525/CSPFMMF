function [Eout,x,dx] = ft(Ein, z, a, b, u, v, lambda, nx, ny)

%This function finds the Fourier Transform by evaluating the integral
%instead of using FFT.  There is propogation.  This is a 2d function.

% a, b = size of first plane (aperture) in meters
% z = distance of propagation in meters
% lambda = wavelength
% nxi, neta = number of points in the first plane (odd)
% u, v = size of plane 2 in meters
% nx, ny = number of points in the second plane (odd)
% Ein = electric field at at first plane

[neta, nxi] = size(Ein);

% dxi = a/(nxi-1);
% deta = b/(neta-1);
% dx = u/(nx-1);
% dy = v/(ny-1);

dxi = a/nxi;
deta = b/neta;
dx = u/nx;
dy = v/ny;

xi = [-(nxi-1)/2:(nxi-1)/2]*dxi;
eta = [-(neta-1)/2:(neta-1)/2]*deta;

x = (-nx/2+0.5:nx/2-0.5)*u/nx;
y = (-ny/2+0.5:ny/2-0.5)*v/ny;
x_xi = x'*xi;
y_eta = y'*eta;

Eout=exp(2*1i*pi*z/lambda)/(1i*lambda*z)*exp(-2*pi*1i*y_eta/lambda/z)*Ein*exp(-2*pi*1i*x_xi/lambda/z)'*dxi*deta;