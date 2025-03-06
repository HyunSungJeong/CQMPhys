function yr = KKi2r (xi, yi, varargin)
% < Description >
%
% yr = KKi2r (xi, yi [,option])
%
% Convert the imaginary part of a causal function into the real part, by
% using the Kramers-Kronig relation. If the input is the real part, then 
% yi = -KKi2r(xr,yr) provides the imaginary part.
% To be more specific, this function computes the Cauchy principal value
% (or principal value integral)
%    yr (xr) = (1/pi) * P.V. \int_{-\infty}^{\infty} dx yi(x) / (x-xr) ,
% where P.V. means the principal value and yi(x), yr(xr) are the functions
% corresponding to the discrete data (xi, yi), (xi, yr), respectively.
% In this sense, the input yi does not need to be real, which might be
% convenient feature in dealing with the correlation functions for multiple
% flavors and non-diagonal elements.
%
% < Input >
% xi : [vector] x points of the imaginary part of a causal function.
% yi : [vector, matrix, or multi-dimensional array] y points of the
%       imaginary part(s) of (a) causal function(s). Each column of yi
%       corresponds to different functions, representing the imaginary part
%       of a function as a piecewise linear function connecting the (x,y)
%       pairs specified xi and the column elements. numel(xi) == size(yi,1)
%       should hold.
%
% < Option >
% 'gflag', .. : [number] It decides how to treat the tail of the function
%       outside of the range of xi. 
%       0 : All zeros. The y values are asummed to drop sharply at the
%           narrow intervals [xi(1)-vsn,xi(1)] and [xi(end),xi(end)+vsn],
%           where vsn = 1e-14.
%       1 : 1/x tail such that yi(1)*xi(1)/abs(x) for left,
%           yi(end)*xi(end)/abs(x) for right.   (Default) 
%       2 : 1/x^2 tail.
%
% < Output >
% yr : [vector, matrix, or multi-dimensional array] The real part of (a)
%       causal function(s), on the x points specifed by xi. It has the same
%       structure as the input 'yi'; each column of yr corresponds to the
%       column of yi with the same column index.
%
% Written by S.Lee (2016)
% Rewritten by S.Lee (Oct.27,2017)
% Documentation updated by S.Lee (Dec. 4,2017)
% Updated by S.Lee (Dec. 6,2017)
% Updated by S.Lee (Jul.02,2020): fixed some mistakes on treating tails;
%       implemented a feature for treating multiple columns of yi.
% Updated by S.Lee (Jun.03,2022): Minor fix to treat multi-dimensional
%       inputs.

% try

gflag = 1; % default
vsn = 1e-14; % In case gflag == 0, very narrow interval to define the sharp drop of yi at the edges
% isx0 = false; % if 0 is included in xi and gflag > 0, skip it temporarily and then interpolate the result for x = 0

while ~isempty(varargin)
    switch varargin{1}
        case 'gflag'
            gflag = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: Unknown option.');
    end
end

% % % sanity check
if isvector(xi)
    xi = xi(:);
else
    error('ERR: input frequency ''xi'' is not a vector.');
end
if isvector(yi)
    yi = yi(:);
end
if numel(xi) ~= size(yi,1)
    error('ERR: numel(xi) ~= size(yi,1)');
end

sz = size(yi); % save the size for future reshaping
yi = yi(:,:); % make yi to matrix; reshape later

if any(diff(xi) < 0)
    disptime('WRN: input ''xi'' is not in ascending order.');
    [xi,ids] = sort(xi,'ascend');
    yi = yi(ids,:);
end
if any(diff(xi) == 0)
    error('ERR: Input ''xi'' is not unique.');
end
if ~any(gflag == (0:2))
    error('ERR: gflag show be 0, 1, or 2.');
end
if gflag > 0
   if xi(1) >= 0
       error('ERR: xi(1) >= 0 -> 1/x^n tail is not well defined.');
   elseif xi(end) <= 0
       error('ERR: xi(end) <= 0 -> 1/x^n tail is not well defined.');
   end
end
% % % % 

dxi = diff(xi);
a = (yi(2:end,:) - yi(1:end-1,:))./dxi; % piecewise slope
b = (-xi(1:end-1).*yi(2:end,:)+xi(2:end).*yi(1:end-1,:))./dxi; % piecewise y intercept

if gflag > 0
    ids = find(xi ~= 0); % to be used for adding the tail contributions
end

yr = zeros(size(yi)); % result

% Contribution of a piecewise linear segment y = (a*x+b) over [x1,x2] to
% the point at x0 is given by the integral of (a*x+b)/(x-x0) over [x1,x2]:
%   a*(x2-x1) + (a*x0+b)*log(abs((x2-x0)/(x1-x0)))  --- (1)
% Since a = (y2-y1)/(x2-x1), y1 = a*x1+b, and y2 = a*x2+b, the sum of the
% first terms will be yi(end)-yi(1).

% contribution from the second term of Eq. (1)
logmat = log((xi(2:end) - xi') ./ (xi(1:end-1) - xi'));
for it = (1:numel(xi)-1)
    % Set components with log(0) to zero.
    logmat(it, it:it+1) = 0;
end
for it = (1:numel(xi))
    yr(it, :) = sum((a * xi(it) + b) .* logmat(:, it), 1);
end
yr(2:end-1,:) = yr(2:end-1,:) + yi(2:end-1,:).*log(abs(dxi(2:end)./dxi(1:end-1)));

% contribution from the first term of Eq. (1)
yr = yr + (yi(end,:)-yi(1,:));


switch gflag
    case 0
% There are no contribution from the outside of the x interval, since they
% are zeros. The sharp edges contribute *only* to the boundary. The below
% contributions to yr(1) and yr(end) can be derived by considering the
% limiting case of x0 -> x1^- or x0 -> x2^+ in Eq. (1).

        yr(1,:) = yr(1,:) + yi(1,:)*log(abs(dxi(1)/vsn));
        yr(end,:) = yr(end,:) - yi(end,:)*log(abs(dxi(end)/vsn));
        
    case 1
% contribution from 1/x tail (= y1*x1/x stretching from the point (x1,y1)
% at the edge) to point at x0:
% \int_{x1}^{inf} dx (y1*x1/x) * (1/(x-x0)) = (y1*x1/x0)*log(abs(x1/(x1-x0))) --- (2)
        
        yr(ids(1:end-1),:) = yr(ids(1:end-1),:) + ...
            xi(end)*log(xi(end)./(xi(end)-xi(ids(1:end-1)))).*(yi(end,:)./xi(ids(1:end-1)));
        
% from the left tail: end <-> 1 / takes opposite sign (-1) to the
% contribution to yr(1:end-1) due to opposite integration interval [-inf, x1]

        yr(ids(2:end),:) = yr(ids(2:end),:) - ...
            xi(1)*log(xi(1)./(xi(1)-xi(ids(2:end)))).*(yi(1,:)./xi(ids(2:end)));
        
% at zero frequency
        yr(xi == 0,:) = yr(xi == 0,:) + (yi(end,:)-yi(1,:));
        
% at the edges: the sum of the second term in Eq.(1) and the term in
% Eq.(2). The divergent terms are cancelled out.
        yr(end,:) = yr(end,:) + yi(end,:)*log(xi(end)/dxi(end));
        yr(1,:) = yr(1,:) - yi(1,:)*log(abs(xi(1)/dxi(1))); % opposite sign simliarly as for yr(2:end)
        
    case 2
% contribution from 1/x^2 tail (= y1*x1^2/x^2 stretching from the point
% (x1,y1) at the edge) to point at x0:
% \int_{x1}^{inf} dx (y1*x1^2/x^2) * (1/(x-x0))
%  = \int_{x1}^{inf} dx (y1*x1^2)*[ -1/x0/x^2 - 1/x0^2/x + 1/x0^2/(x-x0) ]
%  = (-y1*x1^2)*[ 1/x0/x1 + (1/x0^2)*log(abs((x1-x0)/x1)) ]     --- (3)

        yr(ids(1:end-1),:) = yr(ids(1:end-1),:) + ...
            yi(end,:).*(  (-xi(end)./xi(ids(1:end-1))) + ((xi(end)./xi(ids(1:end-1))).^2).* ...
            log(xi(end)./(xi(end)-xi(ids(1:end-1))))  );
        
% from the left tail: end <-> 1 / takes opposite sign (-1) to the
% contribution to yr(1:end-1) due to opposite integration interval [-inf, x1]

        yr(ids(2:end),:) = yr(ids(2:end),:) - ...
            yi(1,:).*(  (-xi(1)./xi(ids(2:end))) + ((xi(1)./xi(ids(2:end))).^2).* ...
            log(xi(1)./(xi(1)-xi(ids(2:end))))  );
        
% at zero frequency; nothing to add
        
% at the edges: the sum of the second term in Eq.(1) and the term in
% Eq.(3). The divergent terms are cancelled out.
        yr(end,:) = yr(end,:) - yi(end,:)*(log(abs(dxi(end)/xi(end))) + 1);
        yr(1,:) = yr(1,:) + yi(1,:)*(log(abs(dxi(1)/xi(1))) + 1); % opposite sign simliarly as for yr(2:end)
end

yr = yr/pi; % factor 1/pi due to the definition of KK relation

if ~isequal(size(yr),sz)
    yr = reshape(yr,sz);
end

% catch e
%     disp2(getReport(e))
%     disp2('Let''s Debug!');
%     keyboard;
% end

end