function A = complex_readmatrix(filename)
	row_str = readlines(filename);
	A = [];
	for i=1:length(row_str)
		boop = erase(row_str(i), ["(", ")"]);
		combined_row = cell2mat(textscan(boop, "%f", "Delimiter", ","))';
		real_row = combined_row(1:2:end);
		imag_row = combined_row(2:2:end);
		A = [A; real_row + 1i * imag_row];
	end
end
