function v = blendv(m1, m2, mask)
    v = m1;
    v(mask,:,:) = m2(mask,:,:);
end