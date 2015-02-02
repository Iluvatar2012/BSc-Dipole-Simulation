#!/usr/bin/env ruby

# basic variables
y_amount 	= 25.0;
t 			= 0.0;

# open new file with argument from console
in_f = File.new(ARGV[0]);

# get additional parameters from console
steps 		= ARGV[1];
steplength 	= ARGV[2];

# convert parameters to numerals
steps		= steps.to_i;
steplength	= steplength.to_i;

# read the first line
line = in_f.readline();

# convert the first line to variable N
number = line.scan(/\d+/);
number.collect! &:to_i;
N = number[0];

# iterate the steplength in y direction
y_tick = Math.sqrt(N/2)/y_amount;

# build a new array for writing to 
densA = Array.new(y_amount);
densB = Array.new(y_amount);

# open a new file
out_f = File.new("density_stat", "a");

# adjust the linenumber
in_f.readline();
in_f.readline();
in_f.readline();
in_f.readline();

# iterate over all steps
i = 0;
while i < steps

	# reset the density values 
	for j in 0..(y_amount-1)
		densA[j] = 0;
		densB[j] = 0;
	end

	# iterate over all particles in the step
	j = 0;
	while j < N

		# read the current line, this is an A particle
		line = in_f.readline();

		# convert the current line to numerals
		numbers = line.scan(/\d+.\d+/);
		numbers.collect! &:to_f;

		# iterate the current index using the y value
		index = numbers[2]/y_tick;
		index = index.floor;

		# increment the respective counter
		densA[index] += 1;

		# read the current line, this is a B particle
		line = in_f.readline();

		# convert the current line to numerals
		numbers = line.scan(/\d+.\d+/);
		numbers.collect! &:to_f;

		# iterate the current index using the y value
		index = numbers[2]/y_tick;
		index = index.floor;

		# increment the respective counter
		densB[index] += 1;

		# increment the counter of the loop
		j += 2;
	end

	# iterate over all entries of both arrays
	for j in 0..(y_amount-1)

		# output to file
		out_f.syswrite(t);
		out_f.syswrite(", ");
		out_f.syswrite((j+0.5)*y_tick);
		out_f.syswrite(", ");
		out_f.syswrite(densA[j]);
		out_f.syswrite(", ");
		out_f.syswrite(densB[j]);
		out_f.syswrite("\n");
	end

	# increase the time value
	t += steplength;

	# increment the loop counter
	i += 1;
end