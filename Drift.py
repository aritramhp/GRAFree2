#!/usr/bin/env python

#-----------------------------------------------------------------
# Compute the Drift - the differece between start and end points 
# of the window
#-----------------------------------------------------------------
def ComputeDrift(gfp_wndw):
	#print x_axis[-1], x_axis[0]
	# print(gfp_wndw[-1])
	# print(gfp_wndw[0])
	# exit()
	drift_x = gfp_wndw[-1][0] - gfp_wndw[0][0]
	drift_y = gfp_wndw[-1][1] - gfp_wndw[0][1]
	return [drift_x, drift_y]


def Drift(GFP, block_size):
	# print(GFP)
	dft = []
	# print(len(GFP), block_size)
	for w in range(len(GFP)-block_size+1):
		window =  min((len(GFP) - w), block_size)
		d = ComputeDrift(GFP[w : w+window])		
		dft.append(d)
		
	return dft
