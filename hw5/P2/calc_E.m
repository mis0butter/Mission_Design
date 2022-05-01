function E1 = calc_E(e1, M1)

f =@(E1,M1,e1) E1 - e1*sin(E1) - M1;
fp =@(E1,e1) 1 - e1*cos(E1);

tol=1e-12;
E1 = M1;
while abs(f(E1,M1,e1)) > tol
    E1 = E1 - f(E1,M1,e1) / fp(E1,e1);
end

end