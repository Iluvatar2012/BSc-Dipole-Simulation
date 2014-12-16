#!/usr/bin/env ruby

# open new file with argument from console
in_f = File.new(ARGV[0]);

# get additional parameters from console
m = ARGV([1]);
G = ARGV([2]);

# read the first line
line = in_f.readline();

# convert the first line to variable N
number = line.scan(/\d+/);
number.collect! &:to_i;
N = number[0];

# open a new file
out_f = File.new("final_stat", "a");

# adjust the linenumber
in_f.readline();
in_f.readline();
in_f.readline();
in_f.readline();
in_f.readline();

# skip N lines, this is the first step which shall be ignored
j=0;
while j<N
	in_f.readline();
	j++;
end

# initiate a couple of variables
sum_A = 0;
sum_B = 0;

sq_sum_A = 0;
sq_sum_B = 0;

# iterate over the next N lines
j=0;

while j<N
	# read the next line
	line = in_f.readline();

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

	# increment j
	j += 2;
end

# normalize all sums
sum_A /= N*0.5;
sum_B /= N*0.5;

sq_sum_A /= N*0.5;
sq_sum_B /= N*0.5;

# calculate the standard deviations for both species
std_dev_A = Math.sqrt(sq_sum_A - sum_A*sum_A);
std_dev_B = Math.sqrt(sq_sum_B - sum_B*sum_B);

# calculate the difference of the sums and the resulting error
diff = Math.abs(sum_A - sum_B);
diff_err = 0;

if (sum_A > sum_B)
	diff_err = std_dev_A - std_dev_B;
else
	diff_err = std_dev_B - std_dev_A;

# output to file
out_f.syswrite(m);
out_f.syswrite(", ");
out_f.syswrite(G);
out_f.syswrite(", ");
out_f.syswrite(diff);
out_f.syswrite(", ");
out_f.syswrite(diff_err);
out_f.syswrite("\n");