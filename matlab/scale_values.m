function scale_values(F, indices, scale)
%SCALE_VALUES Scale selected upper-triangular nonzeros in-place.
    F.scaleValues(indices, scale);
end
