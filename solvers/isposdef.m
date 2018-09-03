function bool = isposdef(A)
[~, p] = chol(A);
bool = p ~= 0;