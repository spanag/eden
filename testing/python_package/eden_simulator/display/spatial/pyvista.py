'''
Display helpers for [PyVista]( https://pyvista.org ).
'''
import numpy as np
import pyvista as pv

def get_neuron_mesh(cell_info):
    '''Get a PyVista mesh, colourable per compartment.
    
    Parameters
    ----------
    
    cell_info: dict[str]
        A `dict` containing `'mesh_vertices'`, `'mesh_faces'` and `'mesh_comp_per_face'` describing a neuron's 3-D mesh, as produced by  by `eden_simulator.experimental.explain_cell`.
        
    Returns
    -------

    mesh: pyvista.UnstructuredGrid
        A 3D mesh which can be assigned scalars per compartment.
    '''
    # NB: get the format of all indices right, or else python will crash hard :D
    # convert the faces to the vtk format: no.of verts, verts, repeat
    mesh_vertices, mesh_faces, mesh_comp_per_face = [cell_info[x] for x in ['mesh_vertices', 'mesh_faces', 'mesh_comp_per_face']]
    mesh_faces = np.array(mesh_faces)
    pv_faces = np.c_[np.full(mesh_faces.shape[0],3),mesh_faces].ravel()
    # print(pv_faces)
    
    # get all the faces per compartment, from compartment per face. Hope it's more efficient than insertion in for loop
    # https://stackoverflow.com/questions/5695349/group-list-by-values/
    import itertools
    seq = list(zip(range(len(mesh_comp_per_face)),mesh_comp_per_face))
    seq.sort(key = lambda x : x[1])
    # print(seq)
    comp_to_faces = [[x[0] for x in data] for (key, data) in itertools.groupby(seq, lambda x : x[1])]
    # print(len(comp_to_faces))

    # Now construct the PolyData as a list of polyhedra (face sets)
    # https://docs.pyvista.org/version/stable/examples/00-load/create-polyhedron
    def GetPolyhedronData(comp_faces): return [(1+len(comp_faces)*4), len(comp_faces)] + list(pv_faces.reshape((-1,4))[comp_faces].ravel())
    pedro = sum([GetPolyhedronData(x) for x in comp_to_faces[:]],[])
    celty = [ pv.CellType.POLYHEDRON] *  len(comp_to_faces)
    grid = pv.UnstructuredGrid(pedro, celty, mesh_vertices)
    
    return grid
