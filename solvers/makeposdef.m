function A = makeposdef(A)

if isposdef(A)
        [L, D] = ldl(A);
        D(D < 0) = abs(D(D < 0));
        A = L*D*L';
end