%Create npoints circles of radius scat_size at random locations on the
%Grid of points zz,yy. Utility function for generating 2d test data
function img = make_dots(zz,yy,npoints,scat_size)
	zmax = max(max(zz));
	zmin = min(min(zz));
	ymax = max(max(yy));
	ymin = min(min(yy));
	[r,c] = size(zz);
	img = zeros(r,c);
	for ii = 1:npoints
		centre = [zmax-zmin;ymax-ymin].*rand(2,1)+[zmin;ymin];
		%keyboard();
		img(vec((zz-centre(1)).^2 + (yy-centre(2)).^2 < scat_size^2)) = 1;
	end;
end
