%% Post-processing function
% Removes intermediate densities

function post_process(x,nelx,nely,volfrac)

[Y,I] = sort(x(:), 'descend');
vt = floor(((volfrac - 0.001) * nelx * nely) / (1 - 0.001));
xd(I(1:vt)) = 1;
xd(I(vt+1:end)) = 0.001;
xd = reshape(xd, nely, nelx);

colormap(gray); imagesc(-xd); axis equal; axis tight; axis off;pause(1e-6);
