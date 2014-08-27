#!/usr/bin/env ruby

# basic variables
bias = 5;

# read open incoming file and read all lines into array
lines = IO.readlines(ARGV[0]);

# convert the first line to variable N
number = lines[0].scan(/\d+/);
number.collect! &:to_i;
N = number[0];

# convert the second line to variable steps
number = lines[1].scan(/\d+/);
number.collect! &:to_i;
steps = number[0];

# open a new file
file = File.new("particle_a_stat", "w");

# iterate over all steps
i = 0;
while i<steps

	# sum over all A particles y-value and their squares
	sum = 0;
	sq_sum = 0;

	# iterate the current time
	numbers = lines[bias + i*N].scan(/\d+.\d+/);
	numbers.collect! &:to_f;
	t = numbers[0];

	# iterate over all A particles
	j=0;
	while j<N
		# convert the current line to numerals, add y position of current particle to sum and squared sum
		numbers = lines[bias + i*N + j].scan(/\d+.\d+/);
		numbers.collect! &:to_f;

		sum += numbers[2];
		sq_sum += numbers[2]*numbers[2];

		# jump to the next A particle
		j+=2;
	end

	# normalize both sums
	sum /= N/2;
	sq_sum /= N/2;

	# compute standard deviation
	std_dev = Math.sqrt(sq_sum-sum*sum);

	# output to file
	file.syswrite(t);
	file.syswrite(", ");
	file.syswrite(sum);
	file.syswrite(", ");
	file.syswrite(std_dev);
	file.syswrite("\n");

	# increment to next step
	i+=1
end