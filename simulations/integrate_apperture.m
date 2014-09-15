%Integrate field, defined on the grid of vertices, over the convex hull
%of the vertices.
%Note that for OCT problem, field is not extracted directly from the data:
%need to compute Ui*(U-Ui) first (easy if defined at same vertices).
%Field is a function of vertex and k
function S = integrate_apperture(vertices, field, varargin)
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
	plane_points = vertices(:,plane_directions); %Don't need the z-coord, it's constant
	
	triangles = delaunay(plane_points); %triangular grid in xy plane
	[ntriangles,t] = size(triangles);
	
	%Integrate over field by breaking it into triangular surfaces
	%We do the integrals for every value of k simultaneously.
	%The projection of each trianglular patch onto the xy plane is the 
	%that computed by the delaunay triangulation. The z coord (height) of
	%each vertex of the triangular patch is given by the value of 'field'
	%at each vertex. The volume enclosed between the patch and its projection
	%onto the xy plane is 1/3 * (sum vertex heights) * (base area in xy plane)
	%This can be deduced by elementary geometry. Note that we can compute
	%the vertex sums for every k simultaneously, and we can allow the vertex
	%sums to be complex (provided we ensure our calculation for the base area
	%returns a positive real).
	[nverts, nk] = size(field);
	S = zeros(1,nk);
	for it = 1:ntriangles
		ia = triangles(it,1);
		ib = triangles(it,2);
		ic = triangles(it,3);
		%Projection of triangle points onto xy plane
		D = [plane_points(ia,1),plane_points(ia,2),0];
		E = [plane_points(ib,1),plane_points(ib,2),0];
		F = [plane_points(ic,1),plane_points(ic,2),0];
		twice_base_area = norm(cross(D-E,F-E));
		%Get field values at the vertices for every k and sum them
		%keyboard();
		sum_field = field(ia,:)+field(ib,:)+field(ic,:);
		%Integral over triangular patch is 1/6*(sum_fields)*twice_base_area
		%Add contribution from this patch onto total
		S = S + (twice_base_area/6)*sum_field;
	end;
	
	%Transpose S to return it as a column vector - don't use ', S may be complex
	S = transpose(S);
	
end
