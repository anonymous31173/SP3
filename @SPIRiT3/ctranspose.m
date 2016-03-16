function res = ctranspose(sp3)
sp3.adjoint = xor(sp3.adjoint,1);
res = sp3;