# BKL251
#Implementor: Sabyasachi Karati

The followings are the descriptions of the available directories:

1. left-to-right/hmult: contains left-to-right Montgomery ladder scalar multiplication
for both the fixed-base and variable-base using hybrid field multiplication module.

2. left-to-right/kmult: contains left-to-right Montgomery ladder scalar multiplication
for both the fixed-base and variable-base using Karatsuba field multiplication module. 

3. left-to-right/smult: contains left-to-right Montgomery ladder scalar multiplication
for both the fixed-base and variable-base using schoolbook field multiplication module.

4. right-to-left/kmult: contains right-to-left Montgomery ladder scalar multiplication
for both the fixed-base and variable-base using Karatsuba field multiplication module.

5. left-to-right/kmult+precomp: contains left-to-right Montgomery ladder scalar
multiplication with precomputation table of the multiples of the base point for the
fixed-base using karatsuba field multiplication module.

The directories gp verification and magma verification contain the relevant scripts to
verify consistency of the scalar multiplications among Kummer line and the associated
elliptic curve.
