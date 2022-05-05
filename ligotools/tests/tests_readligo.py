import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

import readligo as rl


#test 1
def test_read_hdf5():
	rl.read_hdf5(myPath + '/../../' + 'data/L-L1_LOSC_4_V2-1126259446-32.hdf5')


# test 2
def test_loaddata():
	strain, time, chan_dict = rl.loaddata('../data/H-H1_LOSC_4_V2-1126259446-32.hdf5', 'H1')


# test 3
def test_dq_channel_to_seglist():
	strain, time, chan_dict = rl.loaddata(myPath + '/../../' + 'data/L-L1_LOSC_4_V2-1126259446-32.hdf5', 'H1')
	segment_list = rl.dq_channel_to_seglist(chan_dict['CBC_CAT3'])
	print('Number of segments with DQflag CBC_CAT3 = ',len(segment_list))
	assert len(segment_list) == 1


# test 4
def test_FileList():
	fl = rl.FileList()