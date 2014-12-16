#!/usr/bin/env ruby

# open new file with argument from console
in_f = File.new(ARGV[0]);

# read the first line
line = in_f.readline();

# convert the first line to variable N
number = line.scan(/\d+/);
number.collect! &:to_i;
N = number[0];

# read the second line
line = in_f.readline();

# convert the second line to variable steps
number = line.scan(/\d+/);
number.collect! &:to_i;
steps = number[0];

# open a new file
out_f = File.new("particle_stat", "w");

# adjust the linenumber
in_f.readline();
in_f.readline();
in_f.readline();
line = in_f.readline();

# iterate over all steps
i = 0;
while i<steps

	# sum over all A particles y-value and their squares
	sum_A = 0;
	sq_sum_A = 0;

	# sum over all A particles y-value and their squares
	sum_B = 0;
	sq_sum_B = 0;

	# iterate the current time
	numbers = line.scan(/\d+.\d+/);
	numbers.collect! &:to_f;
	t = numbers[0];

	# iterate over all A particles
	j=0;
	while j<N
		# convert the current line to numerals, add y position of current particle to sum and squared sum
		numbers = line.scan(/\d+.\d+/);
		numbers.collect! &:to_f;

		sum_A += numbers[2];
		sq_sum_A += numbers[2]*numbers[2];

		# read the next line
		line = in_f.readline();

		# convert the current line to numerals, add y position of current particle to sum and squared sum
		numbers = line.scan(/\d+.\d+/);
		numbers.collect! &:to_f;

		sum_B += numbers[2];
		sq_sum_B += numbers[2]*numbers[2];

		# read the next line
		line = in_f.readline();

		# jump to the next A particle
		j+=2;
	end

	# normalize both A sums
	sum_A /= N/2;
	sq_sum_A /= N/2;

	# normalize both B sums
	sum_B /= N/2;
	sq_sum_B /= N/2;

	# compute standard deviations
	std_dev_A = Math.sqrt(sq_sum_A-sum_A*sum_A);
	std_dev_B = Math.sqrt(sq_sum_B-sum_B*sum_B);


	# output to file
	out_f.syswrite(t);
	out_f.syswrite(", ");
	out_f.syswrite(sum_A);
	out_f.syswrite(", ");
	out_f.syswrite(std_dev_A);
	out_f.syswrite(", ");
	out_f.syswrite(sum_B);
	out_f.syswrite(", ");
	out_f.syswrite(std_dev_B);
	out_f.syswrite("\n");

	# increment to next step
	i+=1
end