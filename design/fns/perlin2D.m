function In = perlin2D(res,subsmpl)
%% Author:
% Vivek Boominathan
% Rice University
% vivekb@rice.edu

%%
% res - [Wx, Wy] resoln in px of output perlin noise
% subsmpl - [sx, sy] samples between grid points
% ceil(res/subsmpl) - number of grid points

smpl = 1./subsmpl;

%% Permutation lookup
wl = 256; % wrap length

% p = load('p_segu.mat'); p = p.p;
p = randperm(wl);

% To avoid index wrapping, double the table length
p = [p,p];

%% Gradients
ngrads = 16;
g2 = [cos(linspace(0,2*pi*(ngrads-1)/ngrads,ngrads));...
    sin(linspace(0,2*pi*(ngrads-1)/ngrads,ngrads))];

%% Fade or blending fn
fade = @(t) 6*t.^5 - 15*t.^4 + 10*t.^3;

%% Compute perlin noise - vectorized implementation
[x,y] = meshgrid(smpl(1)*(1:res(1)),smpl(2)*(1:res(2)));

% Find the nerest grid point
x0 = floor(x);
x1 = floor(x) + 1;
y0 = floor(y);
y1 = floor(y) + 1;

% Distances (2D)
tx0 = x - x0;
tx1 = x - x1;
ty0 = y - y0;
ty1 = y - y1;

% wrap grids at wrap length
x0 = mod(x0,wl) + 1;
x1 = mod(x1,wl) + 1;
y0 = mod(y0,wl) + 1;
y1 = mod(y1,wl) + 1;

% Calculate the 4 corner grad indices
gi00 = mod(p(x0+p(y0)), ngrads) + 1;
gi01 = mod(p(x0+p(y1)), ngrads) + 1;
gi10 = mod(p(x1+p(y0)), ngrads) + 1;
gi11 = mod(p(x1+p(y1)), ngrads) + 1;

% Calculate the 4 corner grads
g00 = cat(3,reshape(g2(1,gi00),size(gi00)), ...
    reshape(g2(2,gi00),size(gi00)));
g01 = cat(3,reshape(g2(1,gi01),size(gi01)), ...
    reshape(g2(2,gi01),size(gi01)));
g10 = cat(3,reshape(g2(1,gi10),size(gi10)), ...
    reshape(g2(2,gi10),size(gi10)));
g11 = cat(3,reshape(g2(1,gi11),size(gi11)), ...
    reshape(g2(2,gi11),size(gi11)));

% Calculate noise contrib from each 4 corners
n00 = dot(g00,cat(3,tx0,ty0),3);
n01 = dot(g01,cat(3,tx0,ty1),3);
n10 = dot(g10,cat(3,tx1,ty0),3);
n11 = dot(g11,cat(3,tx1,ty1),3);

% Weights of each contribs
w00 = fade(1-tx0).*fade(1-ty0);
w01 = fade(1-tx0).*fade(ty0);
w10 = fade(tx0).*fade(1-ty0);
w11 = fade(tx0).*fade(ty0);

% Net noise
In = w00.*n00 + w01.*n01 + w10.*n10 + w11.*n11;

end