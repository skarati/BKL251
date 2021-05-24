dblKL(x,z) = {
	X = Mod(1,2)*b*(x^2+z^2)^2;
	Z = Mod(1,2)*(x*z)^2;
	[X,Z]
}

diffAdd(x1,z1,x2,z2,x,z) = {
	X = Mod(1,2)*z*(x1*x2+z1*z2)^2;
	Z = Mod(1,2)*x*(x1*z2+x2*z1)^2;
	[X,Z]
}

dblEC(x,z) = {
	X = Mod(1,2)*x^4+Mod(1,2)*z^4*b^4;
	Z = Mod(1,2)*x^2*z^2;

	[X,Z]
}

AddEC(x1,y1,x2,y2) = {
	X = Mod(1,2)*(x1*y2+x2*y1) + Mod(1,2)*x1*x2*(x1+x2);
	Z = Mod(1,2)*(x1^2+x2^2);

	[X,Z]
}

K2E(x,z) = {
	X = Mod(1,2)*b*z;
	Z = Mod(1,2)*x;

	[X,Z]
}

E2K(x,z) = {
	X = Mod(1,2)*b*z;
	Z = Mod(1,2)*x;

	[X,Z]
}

AddECT2(x,z) = {
	X = Mod(1,2)*b^2*x*z; 
	Z = Mod(1,2)*x^2;
	
	[X,Z]
}

\\verification of Eq 8;
P = K2E(Mod(1,2)*x,Mod(1,2)*z);
P = AddECT2(P[1],P[2]);
Q = dblEC(P[1],P[2]);
P = dblKL(x,z);
P = K2E(P[1],P[2]);
R = AddECT2(P[1],P[2]);

check1 = Q[1]*R[2] + Q[2]*R[1];
print("Verification Eq 8: ",lift(check1));


\\verification of 1st equation of Eq 9;
P = AddECT2(Mod(1,2)*x,Mod(1,2)*z);
P = E2K(P[1],P[2]);
Q = dblKL(P[1],P[2]);
R = K2E(Q[1],Q[2]);
R = AddECT2(R[1],R[2]);
S = dblEC(Mod(1,2)*x,Mod(1,2)*z);
check2 = S[1]*R[2] + S[2]*R[1];
print("Verification of 1st equation of Eq 9: ",lift(check2));	


\\verification of 2nd equation of Eq 9;
P = AddEC(Mod(1,2)*x1,Mod(1,2)*y1,Mod(1,2)*x2,Mod(1,2)*(x2+y2));
P = AddECT2(P[1],P[2]);
P = E2K(P[1],P[2]);


P1 = AddECT2(x1,Mod(1,2));
P1 = E2K(P1[1],P1[2]);

P2 = AddECT2(x2,Mod(1,2));
P2 = E2K(P2[1],P2[2]);

R = diffAdd(P1[1],P1[2],P2[1],P2[2],P[1],P[2]);
R = K2E(R[1],R[2]);
R = AddECT2(R[1],R[2]);

Q = AddEC(Mod(1,2)*x1,Mod(1,2)*y1,Mod(1,2)*x2,Mod(1,2)*y2);

check3 = R[1]*Q[2] + R[2]*Q[1];
c5 = polcoeff(check3,5,y1);
c0 = check3 + Mod(1,2)*c5*y1^5;
check3 = Mod(1,2)*c5*(x1*y1+x1^3+b^4)^2*y1+ Mod(1,2)*c0;

c5 = polcoeff(check3,5,y2);
c0 = check3 + Mod(1,2)*c5*y2^5;
check3 = Mod(1,2)*c5*(x2*y2+x2^3+b^4)^2*y2+ Mod(1,2)*c0;

c4 = polcoeff(check3,4,y1);
c0 = check3 + Mod(1,2)*c4*y1^4;
check3 = Mod(1,2)*c4*(x1*y1+x1^3+b^4)^2+ Mod(1,2)*c0;

c4 = polcoeff(check3,4,y2);
c0 = check3 + Mod(1,2)*c4*y2^4;
check3 = Mod(1,2)*c4*(x2*y2+x2^3+b^4)^2+ Mod(1,2)*c0;

c3 = polcoeff(check3,3,y1);
c0 = check3 + Mod(1,2)*c3*y1^3;
check3 = Mod(1,2)*c3*(x1*y1+x1^3+b^4)*y1+ Mod(1,2)*c0;

c3 = polcoeff(check3,3,y2);
c0 = check3 + Mod(1,2)*c3*y2^3;
check3 = Mod(1,2)*c3*(x2*y2+x2^3+b^4)*y2+ Mod(1,2)*c0;

c2 = polcoeff(check3,2,y1);
c0 = check3 + Mod(1,2)*c2*y1^2;
check3 = Mod(1,2)*c2*(x1*y1+x1^3+b^4)+ Mod(1,2)*c0;

c2 = polcoeff(check3,2,y2);
c0 = check3 + Mod(1,2)*c2*y2^2;
check3 = Mod(1,2)*c2*(x2*y2+x2^3+b^4)+ Mod(1,2)*c0;
print("Verification of 2nd equation of Eq 9: ",lift(check2));

