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

# read the amount of steps to bypass
curr_step = ARGV[1].to_i;

# open a new file
file = File.new("data.csv", "w");

# iterate over all particles of one timestep
i = 0
while i<N
	# convert the current line to numerals
	numbers = lines[bias + curr_step*N + i].scan(/\d+.\d+/);
	numbers.collect! &:to_f;

	# write the positions and psi values to file
	file.syswrite(numbers[1]);
	file.syswrite(",");
	file.syswrite(numbers[2]);
	file.syswrite(",");
	file.syswrite(numbers[5]);
	file.syswrite(",");
	file.syswrite(numbers[6]);
	file.syswrite("\n");

	# incerement i
	i += 1;
end 