import sys

print sys.argv[1]

f_in = open(sys.argv[1], 'r')
f_out = open(sys.argv[2], 'w')

last_time = 0;
line_num = 0

for line in f_in:
	if (line[0] != '#'):
		current_time = float(line.split("\t")[0])

#		print line_num, current_time, last_time

		if (current_time < last_time):
			print 'deleted', current_time, line_num
		else:
			f_out.write(line)
			last_time = current_time

	line_num += 1

f_in.close()
f_out.close()