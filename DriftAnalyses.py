#!/usr/bin/env python
import math
import numpy
from numpy.linalg import linalg

#----------------------------------
# Compute the distribution of Mean
#----------------------------------
def ComputeMean(x_axis, y_axis):
	xn_array = numpy.array(x_axis)
	yn_array = numpy.array(y_axis)
	xn_mean = []
	yn_mean = []
	
	mean_x = numpy.mean(xn_array)
	mean_y = numpy.mean(yn_array)
	xn_mean.append(mean_x)
	yn_mean.append(mean_y)
	
	return mean_x, mean_y


#---------------------------------
# Compute the Standard Deviation
#---------------------------------
def ComputeSD(x_axis, y_axis):
	xn_array = numpy.array(x_axis)
	yn_array = numpy.array(y_axis)
	xn_std_dev = []
	yn_std_dev = []
	
	std_dev_x = numpy.std(xn_array)
	std_dev_y = numpy.std(yn_array)
	xn_std_dev.append(std_dev_x)
	yn_std_dev.append(std_dev_y)
			
	return std_dev_x, std_dev_y
	
	
#----------------------------------------------------------------
# Compute Eigen values and Eigen vector of the Covariance Matrix
#----------------------------------------------------------------
def ComputeEigenValues(x_axis, y_axis):
	eigen_value = []
	dominant_angle = []
	
	X = numpy.vstack((x_axis, y_axis))
	# Compute the Covariance Matrix
	#print X
	cov_matrix = numpy.cov(X)
	# Compute eigen value 
	eigval, eigvec = linalg.eig(cov_matrix)
	eigval_list = eigval.tolist()
	eigval = sorted(eigval, key=abs, reverse=True)
	eigen_value.extend(eigval)
	# Compute angle of dominant eigen vector with x-axis
	dominant_pos = eigval_list.index(max(eigval_list, key=abs))
	if eigvec[dominant_pos][0] == 0:
		if eigvec[dominant_pos][1] < 0:
			angle = -90
		else:
			angle = 90
	else:
		angle = math.degrees(math.atan(eigvec[dominant_pos][1] / eigvec[dominant_pos][0]))
	dominant_angle.append(angle)
		
	return eigval, angle


#----------------------------------------------------------------
# Different analyses 
#----------------------------------------------------------------
def DriftAnalyses(x_axis, y_axis):
	#Max & min
	x_min = min(x_axis)
	x_max = max(x_axis)
	y_min = min(y_axis)
	y_max = max(y_axis)
	drift_min = [x_min, y_min]
	drift_max = [x_max, y_max]
	
	# Centroid
	xn_mean, yn_mean = ComputeMean(x_axis, y_axis)
	drift_mean = [xn_mean, yn_mean]
	
	# Standard Deviation
	xn_std_dev, yn_std_dev = ComputeSD(x_axis, y_axis)
	drift_std_dev = [xn_std_dev, yn_std_dev]
	
	# Compute Eigen values and angle of dominant eigen vector of the Covariance Matrix
	eigenvalues, angle_eigenvector = ComputeEigenValues(x_axis, y_axis)
		
	# Start point and end point
	if x_max == x_min:
		mu_S_x = 'inf'
	else:
		mu_S_x = (x_axis[0]-x_min)*1.0/(x_max-x_min)
	if y_max == y_min:
		mu_S_y = 'inf'
	else:
		mu_S_y = (y_axis[0]-y_min)*1.0/(y_max-y_min)
	mu_S = [mu_S_x, mu_S_y]
	if x_max == x_min:
		mu_E_x = 'inf'
	else:
		mu_E_x = (x_axis[-1]-x_min)*1.0/(x_max-x_min)
	if y_max == y_min:
		mu_E_y = 'inf'
	else:
		mu_E_y = (y_axis[-1]-y_min)*1.0/(y_max-y_min)
	mu_E = [mu_E_x, mu_E_y]
	
	# Horizontal and vertical motion
	x_temp = sum([abs(x_axis[i]-x_axis[i+1]) for i in range(len(x_axis)-1)])
	if x_temp == 0:
		mu_HMN = 'inf'
	else:
		mu_HMN = (x_max - x_min)*1.0/x_temp
	y_temp = sum([abs(y_axis[i]-y_axis[i+1]) for i in range(len(y_axis)-1)])
	if y_temp == 0:
		mu_VMN = 'inf'
	else:
		mu_VMN = (y_max - y_min)*1.0/y_temp
	
	# Aspect ratio
	if y_max == y_min and x_max == x_min:
		mu_AR = 'inf'
	else:
		mu_AR = (y_max - y_min)*1.0/(math.sqrt((y_max-y_min)**2 + (x_max-x_min)**2))
	
	# Arcedness and straightness
	if y_max == y_min and x_max == x_min:
		arc = 'inf'
	else:
		arc = len(x_axis) - math.sqrt((y_max-y_min)**2 + (x_max-x_min)**2)
		
	return xn_mean, yn_mean, eigenvalues, angle_eigenvector


#----------------------------------------------------------------------------------------------------
def Analysis(drift, NO_OF_FRAGMENTS):
	# Drift
	drift_feature_vector = []
	drift_xn_mean = []
	drift_yn_mean = []
	drift_eigenvalues = []
	drift_angle_eigenvector = []

	seq_len = len(drift)
	frag_size = round(2*math.log(seq_len, 4)+.5)
	start_end = ComputeStartEnd(frag_size, NO_OF_FRAGMENTS, seq_len)
	# print(start_end)

	dft = list(zip(*drift))
	drift_x, drift_y = list(dft[0]), list(dft[1])
	for start, end in start_end:
		# collect x axis data of Drift	
		frag_x_axis = drift_x[start:end]
		# collect y axis data of Drift
		frag_y_axis = drift_y[start:end]
				
		xn_mean, yn_mean, eigenvalues, angle_eigenvector = DriftAnalyses(frag_x_axis, frag_y_axis)
				
		drift_feature_vector.append([xn_mean, yn_mean, eigenvalues[0], eigenvalues[1], angle_eigenvector])
		
	return drift_feature_vector


#---------------------------------------------------------------------------------------------------
# Compute the start and end point of fragments
#---------------------------------------------------------------------------------------------------
def ComputeStartEnd(frag_size, no_of_frag, seq_len):
	start = 0
	end = frag_size
	extra = frag_size * no_of_frag - seq_len
	overlap = round(extra*1.0/(no_of_frag-1))
	
	start_end = []
	start_end.append([start, end])

	for i in range(1,no_of_frag-1):
		start = end-overlap
		end = min(start + frag_size, seq_len-1)
		extra = frag_size*(no_of_frag-i) - (seq_len-start)
		overlap = round(extra*1.0/(no_of_frag-i-1)+0.5)
		
		start_end.append([start, end])

	start = end-overlap
	end = min(start + frag_size, seq_len-1)
	start_end.append([start, end])

	return start_end
