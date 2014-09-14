%Extract data over a surface
%Use extract_data_opts() to pass optional parameters, see that function
%for details.
function [plane_vertices, Ex_xy_k, plane_indices] = extract_data(fn, varargin)
	if nargin > 1
		opts = varargin{1};
	else
		opts = extract_data_opts();
	end;
	
	%Find all vertices in the plane z=a
	zvals = vertices(:,opts.plane);
	plane_indices = find(zvals == opts.a);
	plane_vertices = vertices(plane_indices,:);
	clear('vertices'); %Free memory, we need it!
	clear('zvals');
	
	%Now extract data at the required vertices
	load(fn, 'camplitudes');
	%Extract required component of E/H, at every wavelength
	Ex_xy_k = camplitudes(plane_indices, 3*opts.EH+opts.component, :);
	clear('camplitudes');
	
	%Todo - is there a way to extract the wavenumbers, or do we just have to know them?
	
	%if there's an output file name, save the data for later use
	if (opts.out ~= '')	
		save('-hdf5', opts.out, 'Ex_xy_k','plane_vertices','plane_indices','plane','a','EH','component');
	end;
end
