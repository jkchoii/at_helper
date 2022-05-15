function [range_src2recv, dir_name] = getSampledRange(start_range,end_range, num_samples)

num_start = length(start_range);
num_end = length(end_range);
assert(num_start == num_end, "error: number of range");

range_src2recv = zeros(num_start * num_samples, 1);
dir_name = string(zeros(num_start, 1));
for idx = 1:num_start
    range_src2recv(1+num_samples*(idx-1) : num_samples*idx) = sort(start_range(idx) + (end_range(idx) - start_range(idx)).*rand(1, num_samples)); % meter
    dir_name(idx) = string(['range' num2str(start_range(idx) ./ 1000) '-' num2str(end_range(idx) ./ 1000) 'km']);
end

end

