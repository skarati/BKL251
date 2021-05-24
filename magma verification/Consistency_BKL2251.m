p:=251;
q:=2^p;
F2 := FiniteField(2);
P<t> := PolynomialRing(F2);
f:= P!(t^251+t^7+t^4+t^2+1);
Z := IntegerRing();
F<z> := ext< F2 | f >;

b:=  F!z^13 + z^9 + z^8 + z^7 + z^2 + z + 1;
E := EllipticCurve([F!1, 0, 0, 0, b^4]);
ordE := 3618502788666131106986593281521497120350295568660281842639600273604847968132;
ord :=  904625697166532776746648320380374280087573892165070460659900068401211992033;

Mapping_binKL_2_EC := function(Pt)
	xt := Pt[2]*b/Pt[1];
	ptw := Points(E,xt);
	return ptw[1];
end function;


Mapping_EC_2_binKL := function(Pt)
	xt := b/Pt[1];
	return [xt,1];
end function;


basePtKL := [F!z^3 + z^2,1];
xt := b/basePtKL[1];
ptw := Points(E,xt);
basePtEC := E!ptw[1];
ptOrd2 := E![0,b^2];

//binKL Arithmetic

KL_double := function(pt)
	xt := b*(pt[1]^2+pt[2]^2)^2;
	zt := (pt[1]*pt[2])^2;
	return [xt,zt];
end function;

KL_diffadd := function(pt1,pt2,P)
	xt := P[2]*(pt1[1]*pt2[1]+pt1[2]*pt2[2])^2;
	zt := P[1]*(pt1[1]*pt2[2]+pt1[2]*pt2[1])^2;
	return [xt,zt];
end function;

//Sanity Check
ScalarMultKL := function(P,n)
	n_bin := Intseq(n,2);
	n_bin := Reverse(n_bin);
	len := #n_bin;
	
	S := P;
	R := KL_double(P);
	for i:= 2 to len do
		if n_bin[i] eq 0 then
			R := KL_diffadd(S,R,P);
			S := KL_double(S);
		elif n_bin[i] eq 1 then
				S := KL_diffadd(S,R,P);
				R := KL_double(R);
		end if;
	end for;
	return [S[1]/S[2],1];
end function;

count:=0;
for n:=1 to 100 do
	ptKL1 := ScalarMultKL(basePtKL,n);
	ptKL2EC := Mapping_binKL_2_EC(ptKL1)+ptOrd2;
	
	PT := Mapping_binKL_2_EC(basePtKL) + ptOrd2;
	ptEC := n*PT;
	count := count + (ptKL2EC[1]+ptEC[1]);	
end for;
print count;


basePtEC:=basePtEC+ptOrd2;

count:=0;
for n:=1 to 100 do
	ptEC2 := n*basePtEC;
	ptEC2KL := Mapping_EC_2_binKL(ptEC2+ptOrd2);

	PT := Mapping_EC_2_binKL(basePtEC+ptOrd2);
	ptKL2 := ScalarMultKL(PT,n);
	count := count + (ptEC2KL[1] + ptKL2[1]);
end for;
print count;

