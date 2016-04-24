import os

new_filenames = os.listdir(os.getcwd() + '/output/')

for i in range(0, len(new_filenames)):
	# print '\nChecking: ' + new_filenames[i]
	orig = [line.rstrip('\n') for line in open(os.getcwd() + '/koes/' + new_filenames[i])]
	new =  [line.rstrip('\n') for line in open(os.getcwd() + '/output/' + new_filenames[i])]
	
	if len(orig) != len(new):
		print new_filenames[i] + " wrong number of lines: %d vs %d" % (len(orig),len(new))
		continue

	for j in range(1, len(orig)):
		# Split on the first '|'
		orig_hold = orig[j].split('|', 1)
		new_hold = new[j].split('|', 1)

		# Split that separates frame number and the distances
		orig_frame_and_dist = orig_hold[1].split('|', 1)
		new_frame_and_dist = new_hold[1].split('|', 1)

		# Frame number is being checked
		orig_frame_num = orig_frame_and_dist[0]
		new_frame_num = new_frame_and_dist[0]

		# Check for atom pattern
		orig_atom_pattern = orig_hold[0].split(':')
		new_atom_pattern = new_hold[0].split(':')

		# Compare distances
		orig_dists = orig_frame_and_dist[1].split('|')
		new_dists = new_frame_and_dist[1][:-1].split('|')

		# Check if atom patterns match
		if orig_atom_pattern == new_atom_pattern:
			for d in range(0, len(new_dists)):
				# Pattern matched, but conflicting distances
				if abs(float(orig_dists[d]) - float(new_dists[d])) > 0.02:
					print new_filenames[i] + " Frame #: " + orig_frame_num + " wrong distance, nmrdump: " + orig_dists[d] + " new: " + new_dists[d] + " " + str(new_atom_pattern)
		else:
			print new_filenames[i] + " Frame #: " + orig_frame_num + " wrong atom pattern, nmrdump: " + str(orig_atom_pattern) + " new: " + str(new_atom_pattern)
			print "Orig",
			for d in orig_dists[::5]:
				print d,
			print ''
			print 'New ',
			for d in new_dists[::5]:
				print d,
			print ''
