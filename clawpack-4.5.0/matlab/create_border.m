function border = create_border(sdir, sval, xe,ye,ze,mappedgrid,manifold);

% Internal matlab routine for Clawpack graphics.

% This is MUCH slower (like 30% slower), but we have to do it this way for
% mapped grids or manifolds, since we need

[xe_like, ye_like, ze_like] = get_xyzlike(xe,ye,ze,sdir);

bh = cell(2,1);

% Draw lines in y-like direction
[ve{1:3}] = get_xyz(xe_like, [ye_like(1) ye_like(end)], ze_like,sdir);
bh{1} = create_mesh(sdir,sval,1,ve{1},ve{2}, ve{3},mappedgrid,manifold);

% Now draw lines in z-like direction
[ve{1:3}] = get_xyz(xe_like, ye_like, [ze_like(1), ze_like(end)],sdir);
bh{2} = create_mesh(sdir,sval,2,ve{1},ve{2},ve{3},mappedgrid,manifold);

border = [bh{1}.xlines bh{2}.ylines];
set(border,'Tag','on');
