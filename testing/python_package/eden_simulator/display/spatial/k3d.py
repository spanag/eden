'''
Straightforward display of neurons with [K3D]( https://k3d-jupyter.org ).
'''
import numpy as np
import scipy
import k3d

from ..animation import GetAutoFps
from . import get_mesh_info, get_verts_faces_per_comp

import pkg_resources

def DeduceColormap(colormap):
	'''Try to deduce a K3D colormap from the passed argument.
	
	See Also
	--------

	:external+k3d:doc:`reference/colormaps`
	
	:py:func:`k3d.helpers.map_colors`
	
	'''
	from k3d import colormaps as k3dmaps
	if isinstance(colormap, str) and len(colormap) > 0:
		k3d_cmap_name = colormap[0].upper()+colormap[1:]
		# fun fact: k3d includes all the colormaps from matplotlib, see https://github.com/K3D-tools/K3D-jupyter/blob/v2.16.1/k3d/colormaps/generate_matplotlib_color_maps.py
		if hasattr(k3dmaps.matplotlib_color_maps, k3d_cmap_name):
			cmap = getattr(k3dmaps.matplotlib_color_maps, k3d_cmap_name)
		# fun fact: k3d includes all the colormaps from matplotlib, see https://github.com/K3D-tools/K3D-jupyter/blob/v2.16.1/k3d/colormaps/generate_paraview_color_maps.py
		elif hasattr(k3dmaps.paraview_color_maps, k3d_cmap_name):
			cmap = getattr(k3dmaps.paraview_color_maps, k3d_cmap_name)
		else:
			raise ValueError("Unknown colormap name `%s`.\nRefer to https://k3d-jupyter.org/reference/colormaps.html for allowed colormap names, or provide your own in the form of list[(0 to 1 index, r, g, b)]" % (colormap))
		return cmap
	else:
		return colormap


def RgbToInt(rgb):
	'''Convert ... x 3 *rgb* array to K3D style packed ints.'''
	# print(25,rgb)
	rgb = np.array(rgb)
	if rgb.shape == (0,): return rgb
	assert( rgb.shape[-1] == 3 )
	# print(np.array2string(rgb, formatter={'float':lambda x:f'{x.hex()}'}))
	rgb =  np.floor(rgb*256).clip(0,256-1).astype('uint32')
	rgb = (rgb[...,0]*256*256) + (rgb[...,1]*256) + (rgb[...,2])
	# print(np.array2string(rgb, formatter={'int': "{0:x}".format}))
	if len(rgb.shape) == 0 : return int(rgb) # some k3d widgets like scalars better as builtin ints
	return rgb

def IntToRgb(hex):
	'''K3D style packed integer array *hex* to ... x 3 *rgb* array.'''
	hex = np.array(hex)
	# print(np.array2string(hex, formatter={'int': "{0:x}".format}))
	assert(np.issubdtype(hex.dtype, int))
	rgb = np.concatenate([ (np.mod(x//pow,256)+.5)/256 for x in [hex[...,None]] for pow in [65536, 256, 1] ], axis=-1)
	# print(np.array2string(rgb, formatter={'float':lambda x:f'{x.hex()}'}))
	return rgb.astype('float32')
	
def plot_neuron(cell_info, comp_values=(.5,.5,.5), time_axis_sec=None, color_map=None, compress_cells=True, **kwargs):
	# print(comp_values.shape)
	'''
	Plot a neuron mesh with K3D, optionally coloured and animated.
	For use with `Plot`.

	Parameters
	---

	vertices: Nx3 float array
		The .
	faces: Mx3 int array
		The .
	comp_per_face: M float vector
		The compartment that each face belongs to.
	comp_values, ndarray of float  K x C or T x K or  T x K x C or K or C (or T x c?) nah
		The values per compartment:
		
		* K x 3:
		
	colormap: 
		. Cannot be used with RGB values.
	
	time_axis_sec: T float array
	need time per sec as well...
	kwargs: dict
		Additional parameters to pass through to `k3d.mesh`.
	
	Returns
	---
	
	'''
	if time_axis_sec is not None and list(time_axis_sec) != list(sorted(time_axis_sec)): raise ValueError('time_axis_sec must be in increasing order')
	
	mesh_vertices, mesh_faces, mesh_comp_per_face, n_comps = get_mesh_info(cell_info)
	mapped_verts, mapped_faces, verts_per_comp, faces_per_comp = get_verts_faces_per_comp(mesh_vertices, mesh_faces, mesh_comp_per_face, n_comps)
	# two modes: mapped, or direct rgb...
	# get colormap from k3d or matplotlib.
	comp_values= np.array(comp_values)
	vs = comp_values.shape
	T = None
	if time_axis_sec is not None:
		time_axis_sec = np.array(time_axis_sec)
		if len(time_axis_sec.shape) != 1: raise ValueError('time_axis_sec should be a one-dimensional vector')
		T = len(time_axis_sec)

	# Check dimensions of value matrix and reshape to time x comps x colorchans
	def mustBeComps(axis):
		if vs[axis] != n_comps: raise ValueError(f'There are {n_comps} compartments, but comp_values.shape[{axis}] == {vs[axis]}')
	def mustBeColor(axis):
		if vs[axis] not in (1, 3): raise ValueError(f'comp_values.shape[{axis}] == {vs[axis]}, should be 1 (for pseudo color) or 3 (for RGB)')
	def mustBeTimee(axis):
		if T is not None and vs[axis] != T:  raise ValueError(f'There are {T} time points defined on time_axis_sec, but comp_values.shape[{axis}] == {vs[axis]}')
	if len(vs) == 3:
		mustBeTimee(0)
		mustBeComps(1)
		mustBeColor(2)
	elif len(vs) == 2:
		if vs[1] == n_comps:
			mustBeTimee(0)
			comp_values = comp_values[:,:,None]
		else:
			mustBeComps(0)
			mustBeColor(1)
			comp_values = comp_values[None,:,:]
	elif len(vs) == 1:
		if vs[0] == n_comps:
			comp_values = comp_values[None,:,None]
		else:
			mustBeColor(0)
			comp_values = comp_values[None,:].repeat(n_comps, axis=0)[None,:,:]
	else:
		raise ValueError(f'Unknown dimensions of comp_values: {vs} when the number of compartments is {n_comps}')

	if time_axis_sec is None and comp_values.shape[0] > 1:
		# perhaps complain instead LATER?
		time_axis_sec = np.arange(comp_values.shape[0])/25.0 # default fps
	elif time_axis_sec is not None:
		if T is not None and comp_values.shape[0] != T:  raise ValueError(f'There are {T} time points defined on time_axis_sec, but {comp_values.shape[0]} time points provided')
	T = comp_values.shape[0]

	suggested_fps = None
	if time_axis_sec is not None and len(time_axis_sec) > 1:
		suggested_fps = GetAutoFps(time_axis_sec)

	if color_map is not None:
		color_map = DeduceColormap(color_map)

	# color_alpha = None not until 'opacities' is supported LATER or so
	if comp_values.shape[2] > 1:
		is_rgb = True
		if color_map is not None:
			raise ValueError('Explicit RGB color was provided along with a colormap, use one or the other')
		# if comp_values.shape[2] > 3:
		#     color_alpha = comp_values[..., -1]
		#     comp_values = comp_values[...,:-1]
		#     assert(len(color_alpha.shape) == 2)
		
	else: is_rgb = False

	# Build the k3d dictionary.
	# When animated, each animated property is a dict of str(time) to np arrays; otherwise it is the np array itself.
	k3d_my_attrs = {}
	# k3d_my_cellattrs = {}
	
	def SetFramed(attrname, values, real_time=None):
		# print(attrname)
		'''Set a frame under steady or animated mode.'''
		if real_time is None:
			k3d_my_attrs[attrname] = values
		else:
			if attrname not in k3d_my_attrs: k3d_my_attrs[attrname] = {}
			k3d_my_attrs[attrname][str(real_time)] = values 
	def SetFrame(frame = 0, real_time = None, compress_cells = False):
		'''Set the k3d value/rgb attributes for an animated or steady frame.'''
		framdata = comp_values[frame,:,:]
		if compress_cells:
			if is_rgb : SetFramed('cell_colors', RgbToInt(framdata).astype(np.uint32 ).flatten().tolist(), real_time)
			else      : SetFramed('cell_attribute',      (framdata).astype(np.float32).flatten().tolist(), real_time)
			if frame == 0 or frame == comp_values.shape[0]-1:
				# need to put some verts as well, otherwise the mesh will appear gray at first even after assigning the vert material?
				if is_rgb : SetFramed('colors'   , [0], real_time)
				else      : SetFramed('attribute', [0], real_time)
		else:
			frammape = mapped_verts(framdata)
			if is_rgb : SetFramed('colors', (RgbToInt(frammape).astype(np.uint32 )), real_time)
			else      : SetFramed('attribute',       (frammape).astype(np.float32), real_time)
			
		# if color_alpha is not None:
		#     SetFramed('opacities',  mapped_verts(color_alpha[frame,:  ] ).astype(np.float32), real_time)
	
	if time_axis_sec is not None:
		for frame, real_time in enumerate(time_axis_sec):
			SetFrame(frame, real_time, compress_cells)
	else:
		SetFrame(0) # don't bother with compressing by cells if it isn't animated
	
	if compress_cells:
		# https://stackoverflow.com/questions/62614152/converting-scipy-sparse-csr-csr-matrix-to-a-list-of-lists
		k3d_my_attrs['custom_data'] = {
			'cells_per_vert': list(verts_per_comp.T.tolil(copy=False).rows),
		}
		for k in list(k3d_my_attrs.keys()):
			if not k.startswith('cell_'): continue
			# print(k3d_my_attrs[k])
			k3d_my_attrs['custom_data'][k] = k3d_my_attrs.pop(k)
		# also set color_range if unset
		color_range = kwargs.pop('color_range', [])
		if not color_range and 'cell_attribute' in k3d_my_attrs['custom_data']:
			color_range = k3d.helpers.check_attribute_color_range(k3d_my_attrs['custom_data']['cell_attribute'])
		kwargs['color_range'] = color_range
	plt_mesh = k3d.mesh(mesh_vertices.astype(np.float32), mesh_faces.astype(np.uint32),
		color_map=color_map, **k3d_my_attrs, **kwargs)
	
	if suggested_fps:
		setattr(plt_mesh, 'suggested_fps', suggested_fps)
	
	return plt_mesh
# testing:
# 	attr color
# 	still anim
# 	cell or not
# plot += plot_neuron(cell_info, sampled_voltage[:][...,None].repeat(1,axis=-1), time_axis_sec=anim_axis[:], color_map='blot', compress_cells=1==1 );
# plt_mesh = None
# plt_mesh = plot_neuron(cell_info)
# plt_mesh = plot_neuron(cell_info, vals);
# plt_mesh = plot_neuron(cell_info, vals[...,None].repeat(3,axis=-1) * [0.9, 2, 0.4]);
# plt_mesh = plot_neuron(cell_info, vals[None,...].repeat(2,axis=0) * np.array([1,0.5])[:,None], time_axis_sec = [0,1]);
# plt_mesh = plot_neuron(cell_info, vals[None,...,None].repeat(3,axis=-1).repeat(2,axis=0) * [2, 0.4, 0.9] * np.array([1,0.5])[:,None,None],time_axis_sec = [0,1]);

def decompress_cells(k3d_mesh):
	# print(k3d_mesh)
	'''Revert per-cell compression to normal mode.'''
	if not getattr(k3d_mesh, 'custom_data', None): return
	custom_data = k3d_mesh.custom_data
	# print(custom_data)
	if 'cells_per_vert' not in custom_data: return
	cells_per_vert = custom_data.pop('cells_per_vert')
	# XXX it's not the full story, we should consult the cell vals instead.
	vvv = [1 for l in cells_per_vert for c in l]
	jjj = [c for l in cells_per_vert for c in l]
	iii = [v for v,l in enumerate(cells_per_vert) for c in l]
	n_cells = max((c for l in cells_per_vert for c in l ))+1
	# print(234, n_verts)
	cpvm = scipy.sparse.coo_matrix((vvv, (jjj,iii)), shape=(n_cells,len(cells_per_vert))) # need to flip dims because numpy thinks that axis=0 is columns
	cpvm /= (np.ones(n_cells)*cpvm).T
	def mapped_verts(v):
		v = np.array(v)
		x =  (v.T @ cpvm).T
		return x
	
	if 'cell_attribute' in custom_data:
		k3d_mesh.attribute = {
			k: mapped_verts(v).astype(np.float32)
			for k,v in custom_data['cell_attribute'].items()
		}
		if not k3d_mesh.color_range:
			k3d_mesh.color_range = k3d.helpers.check_attribute_color_range(k3d_mesh.attribute)
	if 'cell_colors' in custom_data:
		k3d_mesh.colors = {
			k: RgbToInt(mapped_verts(IntToRgb(v))).astype(np.uint32)
			for k,v in custom_data['cell_colors'].items()
		}

enhancement_js_base = pkg_resources.resource_string(__name__, 'k3d_enhancements.js').decode("utf-8")

# set_timebar_inline = '''target.appendChild(timebar_container);'''
set_timebar_absolu = '''{
	let x = timebar_container.style;
	x.position='absolute'; x.width='100%'; x.bottom='0';
	x.zIndex = 42;
	target.appendChild(timebar_container);
	target.querySelectorAll('svg').forEach((x) => {x.style.bottom='40px'});
}''' # TODO set size to calc(100% - 1.2em) to keep them separated perhaps...
class Plot(k3d.Plot):
	'''
	A subclass of `k3d.Plot` with augmented display capabilities.
	
	To be used with extended elements such as `plot_neuron`, along with the usual `k3d` display elements.
	'''
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
	def __iadd__(self, objs):
		"""Add Drawable to plot."""
		for obj in objs:
			# print(obj)
			suggested_fps = getattr(obj, 'suggested_fps', None)
			if suggested_fps:
				self.fps = max(self.fps, suggested_fps)
		return super().__iadd__(obj)
	
	def display(self, **kwargs):
		#"""Show plot inside ipywidgets.Output()."""
		# The way that Jupyter extensions are made now, 
		# the widgets are totally untouchable from JS AFAIK
		# (unless a kind soul widget-side adds a reference from the DOM target to the encloseds js widget, that is)
		# Therefore augmentations such as per-cell compression have to be rolled back for display. 
		#     (and also interactive manipulation of attributes from python side)
			
		# output = widgets.Output()
		# with output:
		for ob in self.objects:
			decompress_cells(ob)
		super().display(**kwargs)
		
	def get_snapshot(self, compression_level=9, voxel_chunks=[], additional_js_code='', **kwargs):

		enhancement_js = enhancement_js_base
		if self.snapshot_type == 'inline':
			enhancement_js = enhancement_js.replace('//[set_timebar]', set_timebar_absolu)
		else: 
			enhancement_js = enhancement_js.replace('//[set_timebar]', set_timebar_absolu)
		
		addditional_js_code = (
			enhancement_js+'''
			AddOuterControls(K3DInstance);
			AddCellPainting(K3DInstance);
			// move the script tag to before (or lose it) because some css of jupyter adds a scroll bar, TODO tell upstream
			// since it's a callback, document.currentScript is null ...
			let target = K3DInstance.getWorld().targetDOMNode;
			target.parentNode.querySelectorAll('script').forEach((x)=>{target.parentNode.insertBefore(x, target);}) // don't prepend or they'll be reversed
			
			'''+additional_js_code
		)
		return super().get_snapshot(compression_level, voxel_chunks, addditional_js_code, **kwargs)
	
	# Additional shortcut for displaying a snapshot
	def show_html(self, snapshot_type='inline', **kwargs):
		'''Show a snapshot of the plot as an ipynb displayable.
		
		Parameters
		---
		
		snapshot_type: str
		
			The snapshot_type to use in capturing.
		
		kwargs: dict
			
			Other arguments to pass to ``get_snapshot()``.
		
		'''
		old_snapshot_type = self.snapshot_type
		self.snapshot_type = snapshot_type
		
		from IPython.display import display, HTML, IFrame
		show = display(HTML(self.get_snapshot(**kwargs)))
		
		self.snapshot_type = old_snapshot_type
		return show

# Ipywidgets plots are nice and all but they are widgets;
# hence if they are just instantiated, whether displayed or not, they are being included in the exported ipynb or html. 
# Their contents are being saved in a much less efficient way than if the snapshot was made and displayed.
# Because we have to instantiate them, we can't get rid of them
# (I tried but something keeps tracking and adding them at the ipykernel comm level)
# Hence: Minimize a plot's content to the least workable unit. 
# Applicable ONLY when the plot is not actually displayed (eg its snapshot is shown instead)
def MinimizePlot(plot, more_objects=[]):
    '''
    Minimise the state of a K3D plot widget.
    
    Useful for publishing notebooks with K3D plots, when a snapshot is preferable to the fragile ipywidgets machinery and uncompressed widget_state.
    TODO
    Parameters
    ---

    plot: k3d.Plot
        The plot to minimise.
    more_objects: list[k3d.]
    '''
    # If sth was removed, put it in more_objects for minimization
    import k3d.objects as obs
    import numpy as np
    from traitlets import TraitError
    type_to_attrs = {
        'Mesh':['color_map','opacity_function','vertices','indices','normals','uvs','texture','attribute','triangles_attribute','volume','colors','custom_data'],
        'Points':['color_map','opacity_function','positions','colors','point_sizes','attribute'],
        'Text2d':['text','position','size']
    }
    for ob in plot.objects+more_objects:
        typ = ob.type
        if typ not in type_to_attrs:
            print(f'Type not supported: {typ}');break
        for attrname in type_to_attrs[typ]:
            attr = getattr(ob, attrname)
            if not isinstance(attr, np.ndarray) and not attr : continue
            if isinstance(attr,dict): fakev = {k:[] for k,v in attr.items()}
            else: fakev = []
            # print(attrname, fakev)
            altvals = [0, None]
            for v in [fakev] + altvals:
                try: 
                    setattr(ob, attrname,v)
                    if v is dict: setattr(ob, attrname,[])
                except TraitError: continue
                break
