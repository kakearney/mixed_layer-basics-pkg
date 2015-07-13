function [U V] = adv_uv(dt,cor,pgx,pgy,U1,V1)
%ADV_UV Adjusts advection values
%
% [U, V] = adv_uv(dt,cor,pgx,pgy,U1,V1);
%
% This function adjusts the post-mixing values of the water advection terms
% U and V based on the Coriolis force and pressure gradients. 
%
% Input variables:
%
%   dt:     time step (s)
%
%   cor:    Coriolis term
%
%   pgx:    pressure gradient in U direction
%
%   pgy:    pressure gradient in V direction
%
%   U1:     advection in the east-west direction (m/s)
%
%   V1:     advection in the north-south direction (m/s)
%
% Output variables:
%
%   U:      adjusted east-west advection (m/s)
%
%   V:      adjusted north-south advection (m/s)

eps = 1.e-7;

k1u = pgx + cor*V1 - eps*U1;
k1v = pgy - cor*U1 - eps*V1;
U2 = U1 + k1u*(dt/2);
V2 = V1 + k1v*(dt/2);
k2u = pgx + cor*V2 - eps*U2;
k2v = pgy -cor*U2 - eps*V2;
U3 = U1 + k2u*(dt/2);
V3 = V1 + k2v*(dt/2);
k3u = pgx + cor*V3 - eps*U3;
k3v = pgy - cor*U3 - eps*V3;
U4 = U1 + k3u*dt;
V4 = V1 + k3v*dt;
k4u = pgx + cor*V4 - eps*V4 ;
k4v = pgy - cor*U4 - eps*U4;

ku = (k1u + 2*k2u + 2*k3u + k4u)/6;
kv = (k1v + 2*k2v + 2*k3v + k4v)/6;

U = U1 + ku*dt;
V = V1 + kv*dt;
