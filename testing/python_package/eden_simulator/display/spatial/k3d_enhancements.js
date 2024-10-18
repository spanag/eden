// Display enhancements for eden_simulator.display.spatial.k3d
function MaybeAddOuterTimeControls(K3DInstance, controls_container){
	
	let time_from, time_to, has_time;
	function GetGuiFromTo(K3DInstance){
		let lil_time_slider = K3DInstance.GUI.controls.controllersMap.time;
		return [lil_time_slider._min, lil_time_slider._max];
	}
	let old_has_gui = K3DInstance.parameters.menuVisibility;
	if(!old_has_gui){
		// Need to temporarily set up the gui just to get time info...
		K3DInstance.setMenuVisibility(true);
		
		let old_canim = K3DInstance.parameters.cameraAnimation; // LATER what if this was intended false?
		K3DInstance.setCameraAnimation(true);
		[time_from, time_to] = GetGuiFromTo(K3DInstance);
		K3DInstance.setCameraAnimation(old_canim);
		
		K3DInstance.setMenuVisibility(old_has_gui);
	}
	else{
		// the gui is already there
		[time_from, time_to] = GetGuiFromTo(K3DInstance);
	}
	if( (time_to - time_from) < 1e-6 ) return;
	
	var play_button = document.createElement('button');
	play_button.innerHTML = '';
	play_button.style='font-size:x-large; line-height:100%;'
	function showPlaying(){
		play_button.innerHTML = '&#x23F8'; // ⏸
	}
	function showPaused(){
		play_button.innerHTML = '&#x23F5;'; // ⏵
	}
	
	play_button.addEventListener('click', (event) => {
		if(K3DInstance.autoPlayed){
			K3DInstance.stopAutoPlay();
			showPaused();
		}
		else{
			K3DInstance.startAutoPlay();
			showPlaying();
		}
	});
	controls_container.appendChild(play_button);

	var time_slider = document.createElement('input');
	time_slider.type='range';
	time_slider.value = time_from;
	time_slider.min = time_from;
	time_slider.max = time_to  ;
	time_slider.step = 'any';
	time_slider.addEventListener('input', (event) => {
		K3DInstance.setTime(event.target.value);
		// for a better scrubbing experience, break the anim loop that has its own internal state. LATER ask upstream to give a way to control that internal state.
		if(K3DInstance.autoPlayed){
			K3DInstance.stopAutoPlay();
			setTimeout(() => {K3DInstance.startAutoPlay();}, 10);
			// and could also check LATER if the user has clicked on `pause` in the interim and ...
		}
	});
	time_slider.style='width:100%;'

	controls_container.appendChild(time_slider);// insertBefore
	
	// trap K3D.setTime to keep track on this control
	let OldSetTime = K3DInstance.setTime;
	function MySetTime(time){
		OldSetTime(time);
		//console.log(time);
		time_slider.value = time;
	}
	
	K3DInstance.setTime = MySetTime;
	function showAsK3d(){
		if(K3DInstance.autoPlayed) showPlaying();
		else showPaused();
	}
	
	// also trap the lil.gui autoplay button to sync this button's state
	// LATER trap whether the gui is dynamically added, as needed
	if(K3DInstance.GUI.controls){
		let lil_autoplay = K3DInstance.GUI.controls.controllersMap.autoPlay;
		let OldPlayPause = lil_autoplay.object.togglePlay;
		function MyTogglePlay(){
			OldPlayPause();
			showAsK3d();
		}
		lil_autoplay.object.togglePlay = MyTogglePlay;	
	}
	showAsK3d();
	
}
function AddFullscreenControls(K3DInstance, controls_container){
	let world = K3DInstance.getWorld();
	let target = world.targetDOMNode;
	
	let enlarge_button = document.createElement('button');
	enlarge_button.style='font-weight:bold; font-size:x-large; line-height:100%; margin-left: auto;';
	enlarge_button.innerHTML='';// ⤢; 📺
	//enlarge_button.appendChild(button_span);
	enlarge_button.title = 'Maximize';
	function showMaximized(){
		enlarge_button.innerHTML = '&#128471;&#xFE0E;'; // 🗗︎&#128471;&#xFE0E; 🗙︎ 73 &times;
	}
	function showRestored(){
		enlarge_button.innerHTML='&#x2922;';// ⤢; 📺  🗖︎ 70
	}
	enlarge_button.addEventListener('click', (event) => {
		K3DInstance.setFullscreen(!K3DInstance.getFullscreen()); 
		showAsK3d();
	});
	const currentWindow = target.ownerDocument.defaultView
        || target.ownerDocument.parentWindow;
	currentWindow.addEventListener('resize', () => { showAsK3d() }); // this runs last as the last added, it's convenient
	controls_container.appendChild(enlarge_button);
	
	function showAsK3d(){
		if(K3DInstance.getFullscreen()) showMaximized();
		else showRestored();
	}
	
	// no need to trap for resize yet, event('resize') handles it
	showAsK3d();
}
function AddOuterControls(K3DInstance){
	//if(K3DInstance.doneEdenCustom) return;
	let world = K3DInstance.getWorld();
	let target = world.targetDOMNode;
	
	let controls_container = document.createElement('span');
	controls_container.style = 'display:flex;';
	{
		let x = controls_container.style;
		x.position = 'absolute';
		x.width = '100%';
		x.bottom = '0';
		x.zIndex = 42;
		target.appendChild(controls_container);
		target.querySelectorAll('svg').forEach( (x) => {
			x.style.bottom = '40px'
		});
	}
	
	MaybeAddOuterTimeControls(K3DInstance, controls_container);
	// LATER if there are no time controls, place the enlarge button at the top right of the target or so
	AddFullscreenControls(K3DInstance, controls_container);

	//K3DInstance.doneEdenCustom = true;
	return;
}

function colorsToFloat32Array(array) {
	const colorsArray = new Float32Array(array.length * 3);

	array.forEach((color, i) => {
		colorsArray[i * 3] = ((color >> 16) & 255) / 255;
		colorsArray[i * 3 + 1] = ((color >> 8) & 255) / 255;
		colorsArray[i * 3 + 2] = (color & 255) / 255;
	});

	return colorsArray;
}
function LerpRgbHex(a,b,f){
	let r1 = (a & 255);
	let r2 = (b & 255);
	let g1 = ((a >> 8) & 255);
	let g2 = ((b >> 8) & 255);
	let b1 = ((a >> 16) & 255);
	let b2 = ((b >> 16) & 255);

	let rf = Math.round(r1 + f * (r2 - r1));
	let gf = Math.round(g1 + f * (g2 - g1));
	let bf = Math.round(b1 + f * (b2 - b1));

	return( (bf << 16) | (gf << 8) | rf );
}
function interpolate(a, b, f, property) {
	let i;
	let interpolated;
	let minLength;
	let maxLength;

	if (property === 'model_matrix') {
		// Note: we don't necessarily have to pass k3d's own THREE if it's just math calculations.
		const matrix = new THREE.Matrix4();
		const translationA = new THREE.Vector3();
		const rotationA = new THREE.Quaternion();
		const scaleA = new THREE.Vector3();
		const translationB = new THREE.Vector3();
		const rotationB = new THREE.Quaternion();
		const scaleB = new THREE.Vector3();

		matrix.set.apply(matrix, a.data);
		matrix.decompose(translationA, rotationA, scaleA);
		matrix.set.apply(matrix, b.data);
		matrix.decompose(translationB, rotationB, scaleB);

		translationA.lerp(translationB, f);
		rotationA.slerp(rotationB, f);
		scaleA.lerp(scaleB, f);

		matrix.compose(translationA, rotationA, scaleA);
		const d = matrix.toArray();

		return {
			data: new Float32Array([
				d[0], d[4], d[8], d[12],
				d[1], d[5], d[9], d[13],
				d[2], d[6], d[10], d[14],
				d[3], d[7], d[11], d[15],
			]),
			shape: a.shape,
		};
	}

	if (typeof (a) === 'string') {
		return (f > 0.5) ? b : a;
	}

	if (typeof (a) === 'boolean') {
		return (f > 0.5) ? b : a;
	}

	if (typeof (a) === 'number') {
		return a + f * (b - a);
	}

	function lerp(a,b,interpolated,f){
		//console.log('lerping', a, b)
		minLength = Math.min(a.length, b.length);
		maxLength = Math.max(a.length, b.length);
		if(property.endsWith('colors')){
			for (i = 0; i < interpolated.length; i++) {
				interpolated[i] = LerpRgbHex(a[i], b[i], f);
			}
		}
		else{
			for (i = 0; i < interpolated.length; i++) {
				interpolated[i] = a[i] + f * (b[i] - a[i]);
			}
		}
	
		if (minLength !== maxLength) {
			for (i = minLength; i < maxLength; i++) {
				interpolated[i] = a[i] || b[i];
			}
		}
	}
	if (a.data) {
		minLength = Math.min(a.data.length, b.data.length);
		maxLength = Math.max(a.data.length, b.data.length);
		interpolated = new a.data.constructor(maxLength);
		lerp(a.data, b.data, interpolated, f);
		return {
			data: interpolated,
			shape: a.shape,
		};
	}

	minLength = Math.min(a.length, b.length);
	maxLength = Math.max(a.length, b.length);
	interpolated = Array(maxLength);
	lerp(a, b, interpolated, f);
	return interpolated;
}
function interpFromJson(json, time, property){
	let anim8 = json[property];

	// if not animated, just use it...
	if(!(anim8 && typeof(anim8.timeSeries) !== 'undefined')) return anim8;
	
	let keypoints = Object.keys(anim8).reduce((p, k) => {
					if (!Number.isNaN(parseFloat(k))) {
						p.push({v: parseFloat(k), k});
					}
					return p;
	}, []).sort((q, w) => q.v - w.v);
	//console.log(keypoints);
	if (time <= keypoints[0].v) {
		return anim8[keypoints[0].k];
	} else if (time >= keypoints[keypoints.length - 1].v) {
		return anim8[keypoints[keypoints.length - 1].k];
	} else {
		for (let i = 0; i < keypoints.length; i++) {
			if (Math.abs(keypoints[i].v - time) < 0.001) {
				return anim8[keypoints[i].k];// why clone? that's why TODO
			}
			if (keypoints[i].v > time && i > 0) {
				let a = keypoints[i - 1].v;
				let b = keypoints[i].v;
				let f = (time - a) / (b - a);

				return interpolate(
					anim8[keypoints[i - 1].k],
					anim8[keypoints[i].k],
					f,
					property,
				);
			}
			// else time >= last keypoint, already checked.
		}
	}
}
function MaybeApplyTime(time, mesh, json, THREE){
	let celldata = json.custom_data;
	if(!celldata.cells_per_vert) return;
	let cells_per_vert = celldata.cells_per_vert;
	// Extract sample logic for a and b.
	// See also interpolateTimeSeries https://github.com/K3D-tools/K3D-jupyter/blob/v2.16.1/js/src/core/lib/timeSeries.js#L236
	// this is how meshes are run normally, see 'create' and 'update': https://github.com/K3D-tools/K3D-jupyter/blob/main/js/src/providers/threejs/objects/MeshStandard.js#L93
	//console.log(celldata)
	if(celldata.cell_colors){
		let animval = interpFromJson(celldata, time, 'cell_colors');
		animval = colorsToFloat32Array(animval);
		mesh.material.setValues({
			color: 0xffffff,
			vertexColors: THREE.VertexColors,
		}); // this has to be done every time if verts are missing anyway
		PaintVerts(mesh.geometry, 'color', animval, 3);
		// LATER update json data also?
	}
	if(celldata.cell_attribute){
		let animval = interpFromJson(celldata, time, 'cell_attribute').slice();
		let miny = json.color_range[0], maxy = json.color_range[1];
		for(let i in animval){
			animval[i] = (animval[i] - miny) / (maxy - miny);
		}
		PaintVerts(mesh.geometry, 'uv', animval, 1);
		// LATER update json data also?
	}
	// could also have done K3DInstance.reload(, changes) instead...? but the material might be missing if no verts are passed...
	function PaintVerts(geometry, attrname, cell_vals, coldep){
		//console.log(cell_vals);
		//console.log(geometry);
		// And map to verts.
		let nVerts = cells_per_vert.length;
		
		// LATER use the modern api i guess instead of writing on the array?
		if(!(geometry.attributes[attrname] && geometry.attributes[attrname].array.length == nVerts*coldep)){
			geometry.setAttribute(attrname, new THREE.BufferAttribute(new Float32Array(nVerts*coldep), coldep));
			
		}
		let data = geometry.attributes[attrname].array;
		//console.log(data);
		data.fill(0);
		// LATER optimise for known values?
		for(let v = 0; v < cells_per_vert.length; v++){
			let cs = cells_per_vert[v]; let cl = cs.length;
			for(let j = 0; j < cs.length; j++){
				let c = cs[j];
				let vid = v * coldep;
				for(let k = 0; k < coldep; k++) data[vid+k] += cell_vals[c*coldep+k] / cl;
			}
		}
		//console.log(data);
		geometry.attributes[attrname].needsUpdate = true; // we'll also need to trigger a re-render when not animated.
	}
	
}
function AddCellPainting(K3DInstance){
	let world = K3DInstance.getWorld();
	//console.log(K3DInstance);
	let OldRefreshAfter = K3DInstance.refreshAfterObjectsChange;
	let THREE = K3DInstance.Provider.THREE;
	function MyRefreshAfter(isUpdate, force){
		const meshes = Object.values(K3DInstance.getWorld().ObjectsById).filter( (ob) => {return ob.type == 'Mesh';} )
		for(const mesh of meshes){
			//console.log(mesh);
			const json = world.ObjectsListJson[mesh.K3DIdentifier];
			if(json.custom_data){
				//console.log(mesh);
				MaybeApplyTime(K3DInstance.parameters.time,mesh, json, THREE);
			}
		}
		
		OldRefreshAfter(isUpdate, force); // only needed for the first?
	}
	
	K3DInstance.refreshAfterObjectsChange = MyRefreshAfter;
	K3DInstance.refreshAfterObjectsChange(true); // use this for updating things TODO
}
