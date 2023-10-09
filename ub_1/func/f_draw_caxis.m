function f_draw_caxis(origloc, sz)
% Function draws coordinate system triad using 3D arrows
% origloc = location vector of the origin defaulT [0 0 0]
% sz = arrow relative lenght i.e. percentage of X-axis length default =0.1 

if nargin<1
    origloc=[0 0 0];
end
if nargin<2
    sz=0.10; % percentage of the arrow length with respect to shaft length
end

% drawing coordinate axes
v=axis; ll=v(2)-v(1);
% end points
ex=origloc+[ sz*ll 0 0]; ey=origloc+[ 0 sz*ll 0]; ez=origloc+[ 0 0 sz*ll];
% drawing 3D arrow
f_arrow3D(origloc, ex, [1 0 0]); % X axis
text(ex(1), ex(2), ex(3), 'X', 'Fontsize',13,'HorizontalAlignment','Left','VerticalAlignment','Bottom'  )
f_arrow3D(origloc, ey, [0 1 0]); % Yaxis
text(ey(1), ey(2), ey(3), 'Y', 'Fontsize',13,'HorizontalAlignment','Right','VerticalAlignment','Bottom'  )
f_arrow3D(origloc, ez, [0 0 1]); % Zaxis
text(ez(1), ez(2), ez(3), 'Z', 'Fontsize',13,'HorizontalAlignment','Right','VerticalAlignment','Top'  )

light; lighting gouraud; shading interp;
