function numeric_vector = struct2vect(struct_vector)
% Convert structure to numeric vector
field_names = fieldnames(struct_vector);
numeric_vector = zeros([length(struct_vector) 1]);
for i=1:length(field_names)
    numeric_vector(i) = double(struct_vector.(field_names{i}));
end
end

