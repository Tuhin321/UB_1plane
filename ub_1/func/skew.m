function [uv]=skew(u)
% [uv]=skew(u)
%
% Function calculates skew-symmetric matrix uv from vector u
%   u  = 3x1 or 1x3 vector
%   uv = 3x3 matrix
%
uv=[  0   -u(3)  u(2)
     u(3)   0   -u(1)
    -u(2)  u(1)   0  ];
%
% Written by Asko Rouvinen LTKK/IMVe 2002