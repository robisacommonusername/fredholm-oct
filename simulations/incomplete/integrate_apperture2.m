%Integrate field, defined on the grid of vertices, over the convex hull
%of the vertices.
%Note that for OCT problem, field is not extracted directly from the data:
%need to compute Ui*(U-Ui) first (easy if defined at same vertices).
%Field is a function of vertex and k
function S = integrate_apperture2(vertices, field, varargin)
	if nargin > 2
		plane = varargin{1};
	else
		plane = 3; %plane z=a by default
	end;
	
	%We will in the comments assume that we are integrating over z=0
	%It is possible to integrate over x=0 or y=0 by specifying the plane
	%variable, but in the comments I'll just refer to the z=0 or xy plane
	directions = [1,2,3];
	plane_directions = directions(directions ~= plane);
	x = vertices(:,plane_directions(1));
	y = vertices(:,plane_directions(2));
	xa = min(x);
	xb = max(x);
	ya = min(y);
	yb = max(y);
	
	[nverts, nk] = size(field);
	S = zeros(nk,1);
	for ik = 1:nk
		S(ik) = dbl_integrate_adaptive_interpolate_irregular(x,y,field(:,ik),xa,xb,ya,yb);
	end;
	
end
