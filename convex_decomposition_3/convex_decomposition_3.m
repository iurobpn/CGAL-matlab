%% Y = convex_decomposition_3(vertices,facets[,vertices2,facets2])
%	this function is a mex file adapting convex_decomposition_3 function
%	from CGAL library to matlab
%
%	Inputs:
%		vertices: list of vertices as Nx3 matrix;
%		facets: list of facets with the format:
%			facets = {facet_1,facet_2,...,facet_m}
%			where:
%				facet_i = X indices(line number of vector 
%				in 'vertices' matrix) of the vectors as a
%				Xx1 vector.
%	
%	Output:
%		Y: cell array of format:
%			Y = {Polyhedron_1, Polyhedron_2, ..., Polyhedron_Z}
%			where the all the polyhedron are the decomposed
%			solution provided by the CGAL function. Each Polyhedron
%			has the form:
%			Polyhedron_i = {vertices,facets}
%			where vertices and facets are in the same format as
%			the inputs.
%
