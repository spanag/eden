'''
Utilities for displaying solid neurons in euclidean space.
'''
import numpy as np
import scipy

def get_verts_faces_per_comp(mesh_vertices, mesh_faces, mesh_comp_per_face, n_comps):
	'''
	Convert the ...

	Parameters
	---
	
	Returns
	---
	
	(mapped_verts, mapped_faces, verts_per_comp, faces_per_comp): tuple
	
		* mapped_verts: scipy.spmatrix
			lslsls
	
		
	
	'''
	n_faces = len(mesh_faces)
	# face_comp is a faces-sized vector with the compartment that each face represents.
	# mesh_faces is a faces x 3 array with the three vertices forming each face.
	# Assemble sparse matrices matching each compartment to irs corresponding faces, and vertices.
	# Use the coordinate format (v,(i,j)), to fill each nonzero (i[k],j[k]) with value v[k]
	# First for mapping each comp in sequence, to each in face_comp in sequence. (a binary matrix)
	faces_per_comp = scipy.sparse.coo_matrix( ( np.ones(n_faces) , (mesh_comp_per_face , range(n_faces)) ) , shape=(n_comps,n_faces) ).tocsr()
	# Now for the more involved part, getting all the vertices per compartment.
	iii = np.repeat(range(n_faces),3) # for each of the compartments, for each vertex of the triangle
	jjj = mesh_faces.flatten() # for each face, for each vertex of the triangle
	vvv = np.ones(len(jjj)) # there is one vertex (or more LATER if vertices are being shared)
	verts_per_face = scipy.sparse.coo_matrix((vvv,(iii,jjj)),shape=(len(mesh_comp_per_face),len(mesh_vertices))).tocsr()
	verts_per_face /= np.ones(n_faces) @ verts_per_face # for each vertex, divide by faces touching the same vertex. could also use bool ops LATER?
	verts_per_comp = faces_per_comp @ verts_per_face
	def mapped_verts(comp_cols): 
		return verts_per_comp.T @ comp_cols
	def mapped_faces(comp_cols): return faces_per_comp.T @ comp_cols
	# for some reason, results may vary by an ulp or so and the effect may or may not not happen on similar verts? consider LATER
	# print(np.array2string(..., formatter={'float':lambda x:f'{x.as_integer_ratio()}'}))
	return mapped_verts, mapped_faces, verts_per_comp, faces_per_comp

def get_mesh_info(cell_info):
	'''
	Convert the 'explain cell' format to mappings from compartments to corresponding vertices and faces of the mesh.

	Parameters
	
	Returns
	---
	
	'''
	mesh_comp_per_face = cell_info['mesh_comp_per_face']
	mesh_vertices = cell_info['mesh_vertices']
	mesh_faces = cell_info['mesh_faces']
	mesh_faces = np.array(mesh_faces) # until polygonal faces at least (if ever?)
	n_comps = len(cell_info['comp_parent']); 
	return mesh_vertices, mesh_faces, mesh_comp_per_face, n_comps
