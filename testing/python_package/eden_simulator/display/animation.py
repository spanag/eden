import numpy as np


def GetAutoFps(anim_axis, max_auto_fps = 60):
	'''Get an appropriate update frequency for a time series.
	'''
	anim_dt = np.diff(anim_axis)
	# the last sample might have minuscule delta as per our sampling logic, ignore if possible
	if len(anim_dt) > 1: anim_dt = anim_dt[:-1]
	# not sure what else to do, pick the min not the mean for now
	anim_typ_dt = np.min(anim_dt)
	anim_typ_fps = 1/anim_typ_dt
	# prefer round numbers where it doesn't matter
	anim_typ_fps = round(anim_typ_fps, 5)
	# and put a cap to avoid too much updating for no real benefit
	anim_typ_fps = min(anim_typ_fps, max_auto_fps) # otherwise pick manually
	return anim_typ_fps

# LATER perhaps do linear filtering instead...? or use exact frames only?
def subsample_trajectories(time_axis_sec, data=[], animation_speed=0.0030, animation_frames_per_second=30):
	'''
	Subsample a set of time series for animated real-time display.
	
	Parameters
	---
	time_axis_sec : np.array[float]
		The timebase for the given time series, in seconds.
	data : array
		A `series` by `time` array of values. Length of `time` must be same as `time_axis_sec`.
	animation_speed : float
		The ratio of how fast the time series is played back, relative to `time_axis_sec`.
	animation_frames_per_second : float
		The desired sampling rate for the animation. The time series are subsampled but not interpolated, thus the true frame rate may vary.
	
	Returns
	---
	samples_picked: np.array[int]
		The samples used from the provided sequence `data`.  Useful for resampling other time series similar to `data`. May be less than expected if there are not enough samples to fill the grid.
	anim_axis:  np.array[float]
		Animation-time for each picked sample.
	sampled_time_axis: np.array[float]
		Original timebase for each picked sample. Equals `time_axis_sec[samples_picked]`.
	sampled_data: list[list]
		Sub-sampled `data`. Equals `data[samples_picked]`.

	'''
	# TODO why should data be a list of arrays?
	if list(time_axis_sec) != list(sorted(time_axis_sec)): raise ValueError('time_axis_sec must be in increasing order')
	
	data_duration = (time_axis_sec[-1] - time_axis_sec[0])
	# anim_duration = data_duration / animation_speed
	sampling_dt = animation_speed / animation_frames_per_second
	'''
			 Care for the samples,
		 instead of linspacing things,
				to avoid moiré.
	
					  ~ the conductor
	'''
	ideal_sample_count = data_duration // sampling_dt
	ideal_samples = list(time_axis_sec[0] + sampling_dt * np.arange(ideal_sample_count))+[time_axis_sec[-1]]
	# Discretise to existing samples
	time_border_sec = (time_axis_sec[:-1]+time_axis_sec[1:])/2
	actual_samples = np.searchsorted(time_border_sec, ideal_samples)
	# Note: np.searchsorted will retun within [0...len], thus we can use the result directly.
	# print(anim_duration, sampling_dt, ideal_sample_count, actual_samples)
	# if there are repeated samples don't complain, just provide the unique set.
	actual_samples = np.unique(actual_samples)
	# np.arange() +
	# print(33, [x for x in data])
	sampled_time_axis = time_axis_sec[actual_samples]
	sampled_data = [ x[actual_samples] for x in data ]
	anim_axis = sampled_time_axis / animation_speed
	return actual_samples, anim_axis, sampled_time_axis, sampled_data
