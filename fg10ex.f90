subroutine fg10ex(q,s,r,w,v,p,C)
implicit double precision (q,s,r,w,v,p)
integer, parameter :: dp = kind(1.d0)
real(dp), intent(in) :: q,s,r,w,v,p
real(dp), intent(out) :: C(0:10,28)
q2=q*q
q3=q2*q
q4=q3*q
q5=q4*q
q6=q5*q
q7=q6*q
q8=q7*q
q9=q8*q
q10=q9*q
q11=q10*q
q12=q11*q
q13=q12*q
q14=q13*q
q15=q14*q
sv=s*v
sr=s*r
s2=s*s
sw=s*w
sp=s*p
v2=v*v
rv=v*r
wv=v*w
vp=v*p
r2=r*r
rw=r*w
rp=r*p
w2=w*w
wp=w*p
p2=p*p
sv2=sv*v
srv=sv*r
s2v=sv*s
swv=sv*w
svp=sv*p
sr2=sr*r
s2r=sr*s
srw=sr*w
srp=sr*p
s2w=s2*w
s2p=s2*p
s3=s2*s
sw2=sw*w
swp=sw*p
sp2=sp*p
v3=v2*v
rv2=v2*r
wv2=v2*w
v2p=v2*p
r2v=rv*r
rwv=rv*w
rvp=rv*p
w2v=wv*w
wvp=wv*p
r3=r2*r
r2w=r2*w
r2p=r2*p
rw2=rw*w
rwp=rw*p
w3=w2*w
sv3=sv2*v
srv2=sv2*r
s2v2=sv2*s
swv2=sv2*w
sv2p=sv2*p
sr2v=srv*r
s2rv=srv*s
srwv=srv*w
srvp=srv*p
s2vp=s2v*p
s3v=s2v*s
sw2v=swv*w
sr3=sr2*r
s2r2=sr2*s
sr2w=sr2*w
sr2p=sr2*p
s2rp=s2r*p
s3r=s2r*s
srw2=srw*w
s2wv=s2w*v
s2rw=s2w*r
s3w=s2w*s
s2w2=s2w*w
s2wp=s2w*p
s3p=s2p*s
s4=s3*s
v4=v3*v
rv3=v3*r
wv3=v3*w
r2v2=rv2*r
rwv2=rv2*w
r3v=r2v*r
r2wv=r2v*w
r4=r3*r
r3w=r3*w
sv4=sv3*v
srv3=sv3*r
s2v3=sv3*s
sr2v2=srv2*r
s2rv2=srv2*s
sr3v=sr2v*r
s2r2v=sr2v*s
sr4=sr3*r
s2r3=sr3*s
s2wv2=s2wv*v
s2rwv=s2wv*r
s3wv=s2wv*s
s2r2w=s2rw*r
s3rw=s2rw*s
s3w2=s3w*w
s3vp=s3p*v
s3rp=s3p*r
s4p=s3p*s
C=0
C(0,1)=&
+1._dp
C(1,2)=&
+1._dp
C(1,4)=&
+1._dp
C(2,1)=&
-1._dp*q3
C(2,2)=&
+1._dp*q2*s
C(2,3)=&
+1._dp
C(2,4)=&
+1._dp*q2*s
C(2,5)=&
-3._dp*q2
C(2,11)=&
+2._dp*q*s
C(2,12)=&
+1._dp*q2
C(2,23)=&
+2._dp*q*s
C(2,24)=&
+1._dp*q2
C(3,2)=&
-2._dp*q3&
+1._dp*q2*v&
+1._dp*q2*r
C(3,3)=&
+3._dp*q2*s
C(3,4)=&
-1._dp*q3&
+1._dp*q2*v&
+1._dp*q2*r
C(3,11)=&
-6._dp*q2&
+2._dp*q*v&
+2._dp*q*r
C(3,13)=&
+1._dp*q2
C(3,15)=&
+1._dp*q2
C(3,17)=&
+6._dp*q*s
C(3,18)=&
+3._dp*q2
C(3,23)=&
-3._dp*q2&
+2._dp*q*v&
+2._dp*q*r
C(3,25)=&
+1._dp*q2
C(3,27)=&
+1._dp*q2
C(4,1)=&
+2._dp*q6&
-1._dp*q5*v&
-1._dp*q5*r
C(4,2)=&
-2._dp*q5*s&
+1._dp*q4*sv&
+1._dp*q4*sr&
+3._dp*q2*w
C(4,3)=&
-5._dp*q3&
+4._dp*q2*v&
+4._dp*q2*r&
+3._dp*q4*s2
C(4,4)=&
-2._dp*q5*s&
+1._dp*q4*sv&
+1._dp*q4*sr&
+3._dp*q2*w
C(4,5)=&
+12._dp*q5&
-5._dp*q4*v&
-5._dp*q4*r
C(4,7)=&
-1._dp*q5
C(4,9)=&
-1._dp*q5
C(4,11)=&
-10._dp*q4*s&
+4._dp*q3*sv&
+4._dp*q3*sr&
+6._dp*q*w
C(4,12)=&
-2._dp*q5&
+1._dp*q4*v&
+1._dp*q4*r
C(4,13)=&
+1._dp*q4*s
C(4,14)=&
+3._dp*q2
C(4,15)=&
+1._dp*q4*s
C(4,17)=&
-15._dp*q2&
+8._dp*q*v&
+8._dp*q*r&
+12._dp*q3*s2
C(4,18)=&
+6._dp*q4*s
C(4,19)=&
+4._dp*q2
C(4,21)=&
+4._dp*q2
C(4,23)=&
-10._dp*q4*s&
+4._dp*q3*sv&
+4._dp*q3*sr&
+6._dp*q*w
C(4,24)=&
-2._dp*q5&
+1._dp*q4*v&
+1._dp*q4*r
C(4,25)=&
+1._dp*q4*s
C(4,26)=&
+3._dp*q2
C(4,27)=&
+1._dp*q4*s
C(5,1)=&
-6._dp*q5*w
C(5,2)=&
+4._dp*q6&
+3._dp*q2*p&
-4._dp*q5*v&
+1._dp*q4*v2&
-7._dp*q5*r&
+2._dp*q4*rv&
+1._dp*q4*r2&
+9._dp*q4*sw
C(5,3)=&
-21._dp*q5*s&
+15._dp*q4*sv&
+15._dp*q4*sr&
+15._dp*q2*w
C(5,4)=&
+2._dp*q6&
+3._dp*q2*p&
-3._dp*q5*v&
+1._dp*q4*v2&
-6._dp*q5*r&
+2._dp*q4*rv&
+1._dp*q4*r2&
+9._dp*q4*sw
C(5,5)=&
-30._dp*q4*w
C(5,8)=&
-6._dp*q5
C(5,11)=&
+24._dp*q5&
+6._dp*q*p&
-20._dp*q4*v&
+4._dp*q3*v2&
-35._dp*q4*r&
+8._dp*q3*rv&
+4._dp*q3*r2&
+36._dp*q3*sw
C(5,12)=&
+9._dp*q4*w
C(5,13)=&
-7._dp*q5&
+2._dp*q4*v&
+2._dp*q4*r
C(5,14)=&
+9._dp*q4*s
C(5,15)=&
-4._dp*q5&
+2._dp*q4*v&
+2._dp*q4*r
C(5,16)=&
+3._dp*q2
C(5,17)=&
-105._dp*q4*s&
+60._dp*q3*sv&
+60._dp*q3*sr&
+30._dp*q*w
C(5,18)=&
-21._dp*q5&
+15._dp*q4*v&
+15._dp*q4*r
C(5,19)=&
+15._dp*q4*s
C(5,20)=&
+15._dp*q2
C(5,21)=&
+15._dp*q4*s
C(5,23)=&
+12._dp*q5&
+6._dp*q*p&
-15._dp*q4*v&
+4._dp*q3*v2&
-30._dp*q4*r&
+8._dp*q3*rv&
+4._dp*q3*r2&
+36._dp*q3*sw
C(5,24)=&
+9._dp*q4*w
C(5,25)=&
-6._dp*q5&
+2._dp*q4*v&
+2._dp*q4*r
C(5,26)=&
+9._dp*q4*s
C(5,27)=&
-3._dp*q5&
+2._dp*q4*v&
+2._dp*q4*r
C(5,28)=&
+3._dp*q2
C(6,1)=&
-4._dp*q9&
-9._dp*q5*p&
+4._dp*q8*v&
-1._dp*q7*v2&
+13._dp*q8*r&
-2._dp*q7*rv&
-1._dp*q7*r2&
-9._dp*q7*sw
C(6,2)=&
+4._dp*q8*s&
+18._dp*q4*sp&
-4._dp*q7*sv&
+1._dp*q6*sv2&
-13._dp*q7*sr&
+2._dp*q6*srv&
+1._dp*q6*sr2&
-30._dp*q5*w&
+15._dp*q4*wv&
+15._dp*q4*rw&
+9._dp*q6*s2w
C(6,3)=&
+25._dp*q6&
+18._dp*q2*p&
-40._dp*q5*v&
+16._dp*q4*v2&
-58._dp*q5*r&
+32._dp*q4*rv&
+16._dp*q4*r2&
-30._dp*q7*s2&
+15._dp*q6*s2v&
+15._dp*q6*s2r&
+99._dp*q4*sw
C(6,4)=&
+4._dp*q8*s&
+18._dp*q4*sp&
-4._dp*q7*sv&
+1._dp*q6*sv2&
-13._dp*q7*sr&
+2._dp*q6*srv&
+1._dp*q6*sr2&
-21._dp*q5*w&
+15._dp*q4*wv&
+15._dp*q4*rw&
+9._dp*q6*s2w
C(6,5)=&
-36._dp*q8&
-45._dp*q4*p&
+32._dp*q7*v&
-7._dp*q6*v2&
+104._dp*q7*r&
-14._dp*q6*rv&
-7._dp*q6*r2&
-63._dp*q6*sw
C(6,6)=&
-9._dp*q7*w
C(6,7)=&
+13._dp*q8&
-2._dp*q7*v&
-2._dp*q7*r
C(6,8)=&
-9._dp*q7*s
C(6,9)=&
+4._dp*q8&
-2._dp*q7*v&
-2._dp*q7*r
C(6,10)=&
-9._dp*q5
C(6,11)=&
+32._dp*q7*s&
+72._dp*q3*sp&
-28._dp*q6*sv&
+6._dp*q5*sv2&
-91._dp*q6*sr&
+12._dp*q5*srv&
+6._dp*q5*sr2&
-150._dp*q4*w&
+60._dp*q3*wv&
+60._dp*q3*rw&
+54._dp*q5*s2w
C(6,12)=&
+4._dp*q8&
+18._dp*q4*p&
-4._dp*q7*v&
+1._dp*q6*v2&
-13._dp*q7*r&
+2._dp*q6*rv&
+1._dp*q6*r2&
+18._dp*q6*sw
C(6,13)=&
-13._dp*q7*s&
+2._dp*q6*sv&
+2._dp*q6*sr&
+15._dp*q4*w
C(6,14)=&
-30._dp*q5&
+15._dp*q4*v&
+15._dp*q4*r&
+9._dp*q6*s2
C(6,15)=&
-4._dp*q7*s&
+2._dp*q6*sv&
+2._dp*q6*sr&
+15._dp*q4*w
C(6,16)=&
+18._dp*q4*s
C(6,17)=&
+150._dp*q5&
+36._dp*q*p&
-200._dp*q4*v&
+64._dp*q3*v2&
-290._dp*q4*r&
+128._dp*q3*rv&
+64._dp*q3*r2&
-210._dp*q6*s2&
+90._dp*q5*s2v&
+90._dp*q5*s2r&
+396._dp*q3*sw
C(6,18)=&
-60._dp*q7*s&
+30._dp*q6*sv&
+30._dp*q6*sr&
+99._dp*q4*w
C(6,19)=&
-58._dp*q5&
+32._dp*q4*v&
+32._dp*q4*r&
+15._dp*q6*s2
C(6,20)=&
+99._dp*q4*s
C(6,21)=&
-40._dp*q5&
+32._dp*q4*v&
+32._dp*q4*r&
+15._dp*q6*s2
C(6,22)=&
+18._dp*q2
C(6,23)=&
+32._dp*q7*s&
+72._dp*q3*sp&
-28._dp*q6*sv&
+6._dp*q5*sv2&
-91._dp*q6*sr&
+12._dp*q5*srv&
+6._dp*q5*sr2&
-105._dp*q4*w&
+60._dp*q3*wv&
+60._dp*q3*rw&
+54._dp*q5*s2w
C(6,24)=&
+4._dp*q8&
+18._dp*q4*p&
-4._dp*q7*v&
+1._dp*q6*v2&
-13._dp*q7*r&
+2._dp*q6*rv&
+1._dp*q6*r2&
+18._dp*q6*sw
C(6,25)=&
-13._dp*q7*s&
+2._dp*q6*sv&
+2._dp*q6*sr&
+15._dp*q4*w
C(6,26)=&
-21._dp*q5&
+15._dp*q4*v&
+15._dp*q4*r&
+9._dp*q6*s2
C(6,27)=&
-4._dp*q7*s&
+2._dp*q6*sv&
+2._dp*q6*sr&
+15._dp*q4*w
C(6,28)=&
+18._dp*q4*s
C(7,1)=&
-36._dp*q7*sp&
+60._dp*q8*w&
-30._dp*q7*wv&
-30._dp*q7*rw
C(7,2)=&
-8._dp*q9&
-57._dp*q5*p&
+12._dp*q8*v&
+33._dp*q4*vp&
-6._dp*q7*v2&
+1._dp*q6*v3&
+60._dp*q8*r&
+33._dp*q4*rp&
-36._dp*q7*rv&
+3._dp*q6*rv2&
-30._dp*q7*r2&
+3._dp*q6*r2v&
+1._dp*q6*r3&
+45._dp*q6*s2p&
-108._dp*q7*sw&
+54._dp*q6*swv&
+54._dp*q6*srw&
+45._dp*q4*w2
C(7,3)=&
+144._dp*q8*s&
+189._dp*q4*sp&
-198._dp*q7*sv&
+63._dp*q6*sv2&
-324._dp*q7*sr&
+126._dp*q6*srv&
+63._dp*q6*sr2&
-267._dp*q5*w&
+210._dp*q4*wv&
+210._dp*q4*rw&
+252._dp*q6*s2w
C(7,4)=&
-4._dp*q9&
-39._dp*q5*p&
+8._dp*q8*v&
+33._dp*q4*vp&
-5._dp*q7*v2&
+1._dp*q6*v3&
+38._dp*q8*r&
+33._dp*q4*rp&
-34._dp*q7*rv&
+3._dp*q6*rv2&
-29._dp*q7*r2&
+3._dp*q6*r2v&
+1._dp*q6*r3&
+45._dp*q6*s2p&
-90._dp*q7*sw&
+54._dp*q6*swv&
+54._dp*q6*srw&
+45._dp*q4*w2
C(7,5)=&
-252._dp*q6*sp&
+480._dp*q7*w&
-210._dp*q6*wv&
-210._dp*q6*rw
C(7,6)=&
-36._dp*q7*p
C(7,7)=&
-30._dp*q7*w
C(7,8)=&
+60._dp*q8&
-30._dp*q7*v&
-30._dp*q7*r
C(7,9)=&
-30._dp*q7*w
C(7,10)=&
-36._dp*q7*s
C(7,11)=&
-72._dp*q8&
-285._dp*q4*p&
+96._dp*q7*v&
+132._dp*q3*vp&
-42._dp*q6*v2&
+6._dp*q5*v3&
+480._dp*q7*r&
+132._dp*q3*rp&
-252._dp*q6*rv&
+18._dp*q5*rv2&
-210._dp*q6*r2&
+18._dp*q5*r2v&
+6._dp*q5*r3&
+270._dp*q5*s2p&
-756._dp*q6*sw&
+324._dp*q5*swv&
+324._dp*q5*srw&
+180._dp*q3*w2
C(7,12)=&
+90._dp*q6*sp&
-108._dp*q7*w&
+54._dp*q6*wv&
+54._dp*q6*rw
C(7,13)=&
+60._dp*q8&
+33._dp*q4*p&
-36._dp*q7*v&
+3._dp*q6*v2&
-60._dp*q7*r&
+6._dp*q6*rv&
+3._dp*q6*r2&
+54._dp*q6*sw
C(7,14)=&
-108._dp*q7*s&
+54._dp*q6*sv&
+54._dp*q6*sr&
+90._dp*q4*w
C(7,15)=&
+12._dp*q8&
+33._dp*q4*p&
-12._dp*q7*v&
+3._dp*q6*v2&
-36._dp*q7*r&
+6._dp*q6*rv&
+3._dp*q6*r2&
+54._dp*q6*sw
C(7,16)=&
-57._dp*q5&
+33._dp*q4*v&
+33._dp*q4*r&
+45._dp*q6*s2
C(7,17)=&
+1152._dp*q7*s&
+756._dp*q3*sp&
-1386._dp*q6*sv&
+378._dp*q5*sv2&
-2268._dp*q6*sr&
+756._dp*q5*srv&
+378._dp*q5*sr2&
-1335._dp*q4*w&
+840._dp*q3*wv&
+840._dp*q3*rw&
+1512._dp*q5*s2w
C(7,18)=&
+144._dp*q8&
+189._dp*q4*p&
-198._dp*q7*v&
+63._dp*q6*v2&
-324._dp*q7*r&
+126._dp*q6*rv&
+63._dp*q6*r2&
+504._dp*q6*sw
C(7,19)=&
-324._dp*q7*s&
+126._dp*q6*sv&
+126._dp*q6*sr&
+210._dp*q4*w
C(7,20)=&
-267._dp*q5&
+210._dp*q4*v&
+210._dp*q4*r&
+252._dp*q6*s2
C(7,21)=&
-198._dp*q7*s&
+126._dp*q6*sv&
+126._dp*q6*sr&
+210._dp*q4*w
C(7,22)=&
+189._dp*q4*s
C(7,23)=&
-36._dp*q8&
-195._dp*q4*p&
+64._dp*q7*v&
+132._dp*q3*vp&
-35._dp*q6*v2&
+6._dp*q5*v3&
+304._dp*q7*r&
+132._dp*q3*rp&
-238._dp*q6*rv&
+18._dp*q5*rv2&
-203._dp*q6*r2&
+18._dp*q5*r2v&
+6._dp*q5*r3&
+270._dp*q5*s2p&
-630._dp*q6*sw&
+324._dp*q5*swv&
+324._dp*q5*srw&
+180._dp*q3*w2
C(7,24)=&
+90._dp*q6*sp&
-90._dp*q7*w&
+54._dp*q6*wv&
+54._dp*q6*rw
C(7,25)=&
+38._dp*q8&
+33._dp*q4*p&
-34._dp*q7*v&
+3._dp*q6*v2&
-58._dp*q7*r&
+6._dp*q6*rv&
+3._dp*q6*r2&
+54._dp*q6*sw
C(7,26)=&
-90._dp*q7*s&
+54._dp*q6*sv&
+54._dp*q6*sr&
+90._dp*q4*w
C(7,27)=&
+8._dp*q8&
+33._dp*q4*p&
-10._dp*q7*v&
+3._dp*q6*v2&
-34._dp*q7*r&
+6._dp*q6*rv&
+3._dp*q6*r2&
+54._dp*q6*sw
C(7,28)=&
-39._dp*q5&
+33._dp*q4*v&
+33._dp*q4*r&
+45._dp*q6*s2
C(8,1)=&
+8._dp*q12&
+153._dp*q8*p&
-12._dp*q11*v&
-99._dp*q7*vp&
+6._dp*q10*v2&
-1._dp*q9*v3&
-120._dp*q11*r&
-99._dp*q7*rp&
+66._dp*q10*rv&
-3._dp*q9*rv2&
+60._dp*q10*r2&
-3._dp*q9*r2v&
-1._dp*q9*r3&
-45._dp*q9*s2p&
+108._dp*q10*sw&
-54._dp*q9*swv&
-54._dp*q9*srw&
-135._dp*q7*w2
C(8,2)=&
-8._dp*q11*s&
-414._dp*q7*sp&
+12._dp*q10*sv&
+243._dp*q6*svp&
-6._dp*q9*sv2&
+1._dp*q8*sv3&
+120._dp*q10*sr&
+243._dp*q6*srp&
-66._dp*q9*srv&
+3._dp*q8*srv2&
-60._dp*q9*sr2&
+3._dp*q8*sr2v&
+1._dp*q8*sr3&
+45._dp*q8*s3p&
+252._dp*q8*w&
+189._dp*q4*wp&
-252._dp*q7*wv&
+63._dp*q6*wv2&
-414._dp*q7*rw&
+126._dp*q6*rwv&
+63._dp*q6*r2w&
-108._dp*q9*s2w&
+54._dp*q8*s2wv&
+54._dp*q8*s2rw&
+297._dp*q6*sw2
C(8,3)=&
-152._dp*q9&
-513._dp*q5*p&
+354._dp*q8*v&
+432._dp*q4*vp&
-267._dp*q7*v2&
+64._dp*q6*v3&
+795._dp*q8*r&
+432._dp*q4*rp&
-894._dp*q7*rv&
+192._dp*q6*rv2&
-627._dp*q7*r2&
+192._dp*q6*r2v&
+64._dp*q6*r3&
+252._dp*q10*s2&
+864._dp*q6*s2p&
-252._dp*q9*s2v&
+63._dp*q8*s2v2&
-504._dp*q9*s2r&
+126._dp*q8*s2rv&
+63._dp*q8*s2r2&
-2286._dp*q7*sw&
+1566._dp*q6*swv&
+1566._dp*q6*srw&
+252._dp*q8*s3w&
+675._dp*q4*w2
C(8,4)=&
-8._dp*q11*s&
-324._dp*q7*sp&
+12._dp*q10*sv&
+243._dp*q6*svp&
-6._dp*q9*sv2&
+1._dp*q8*sv3&
+120._dp*q10*sr&
+243._dp*q6*srp&
-66._dp*q9*srv&
+3._dp*q8*srv2&
-60._dp*q9*sr2&
+3._dp*q8*sr2v&
+1._dp*q8*sr3&
+45._dp*q8*s3p&
+144._dp*q8*w&
+189._dp*q4*wp&
-198._dp*q7*wv&
+63._dp*q6*wv2&
-360._dp*q7*rw&
+126._dp*q6*rwv&
+63._dp*q6*r2w&
-108._dp*q9*s2w&
+54._dp*q8*s2wv&
+54._dp*q8*s2rw&
+297._dp*q6*sw2
C(8,5)=&
+96._dp*q11&
+1224._dp*q7*p&
-132._dp*q10*v&
-693._dp*q6*vp&
+60._dp*q9*v2&
-9._dp*q8*v3&
-1320._dp*q10*r&
-693._dp*q6*rp&
+660._dp*q9*rv&
-27._dp*q8*rv2&
+600._dp*q9*r2&
-27._dp*q8*r2v&
-9._dp*q8*r3&
-405._dp*q8*s2p&
+1080._dp*q9*sw&
-486._dp*q8*swv&
-486._dp*q8*srw&
-945._dp*q6*w2
C(8,6)=&
-90._dp*q9*sp&
+108._dp*q10*w&
-54._dp*q9*wv&
-54._dp*q9*rw
C(8,7)=&
-120._dp*q11&
-99._dp*q7*p&
+66._dp*q10*v&
-3._dp*q9*v2&
+120._dp*q10*r&
-6._dp*q9*rv&
-3._dp*q9*r2&
-54._dp*q9*sw
C(8,8)=&
+108._dp*q10*s&
-54._dp*q9*sv&
-54._dp*q9*sr&
-270._dp*q7*w
C(8,9)=&
-12._dp*q11&
-99._dp*q7*p&
+12._dp*q10*v&
-3._dp*q9*v2&
+66._dp*q10*r&
-6._dp*q9*rv&
-3._dp*q9*r2&
-54._dp*q9*sw
C(8,10)=&
+153._dp*q8&
-99._dp*q7*v&
-99._dp*q7*r&
-45._dp*q9*s2
C(8,11)=&
-88._dp*q10*s&
-2898._dp*q6*sp&
+120._dp*q9*sv&
+1458._dp*q5*svp&
-54._dp*q8*sv2&
+8._dp*q7*sv3&
+1200._dp*q9*sr&
+1458._dp*q5*srp&
-594._dp*q8*srv&
+24._dp*q7*srv2&
-540._dp*q8*sr2&
+24._dp*q7*sr2v&
+8._dp*q7*sr3&
+360._dp*q7*s3p&
+2016._dp*q7*w&
+756._dp*q3*wp&
-1764._dp*q6*wv&
+378._dp*q5*wv2&
-2898._dp*q6*rw&
+756._dp*q5*rwv&
+378._dp*q5*r2w&
-972._dp*q8*s2w&
+432._dp*q7*s2wv&
+432._dp*q7*s2rw&
+1782._dp*q5*sw2
C(8,12)=&
-8._dp*q11&
-414._dp*q7*p&
+12._dp*q10*v&
+243._dp*q6*vp&
-6._dp*q9*v2&
+1._dp*q8*v3&
+120._dp*q10*r&
+243._dp*q6*rp&
-66._dp*q9*rv&
+3._dp*q8*rv2&
-60._dp*q9*r2&
+3._dp*q8*r2v&
+1._dp*q8*r3&
+135._dp*q8*s2p&
-216._dp*q9*sw&
+108._dp*q8*swv&
+108._dp*q8*srw&
+297._dp*q6*w2
C(8,13)=&
+120._dp*q10*s&
+243._dp*q6*sp&
-66._dp*q9*sv&
+3._dp*q8*sv2&
-120._dp*q9*sr&
+6._dp*q8*srv&
+3._dp*q8*sr2&
-414._dp*q7*w&
+126._dp*q6*wv&
+126._dp*q6*rw&
+54._dp*q8*s2w
C(8,14)=&
+252._dp*q8&
+189._dp*q4*p&
-252._dp*q7*v&
+63._dp*q6*v2&
-414._dp*q7*r&
+126._dp*q6*rv&
+63._dp*q6*r2&
-108._dp*q9*s2&
+54._dp*q8*s2v&
+54._dp*q8*s2r&
+594._dp*q6*sw
C(8,15)=&
+12._dp*q10*s&
+243._dp*q6*sp&
-12._dp*q9*sv&
+3._dp*q8*sv2&
-66._dp*q9*sr&
+6._dp*q8*srv&
+3._dp*q8*sr2&
-252._dp*q7*w&
+126._dp*q6*wv&
+126._dp*q6*rw&
+54._dp*q8*s2w
C(8,16)=&
-414._dp*q7*s&
+243._dp*q6*sv&
+243._dp*q6*sr&
+45._dp*q8*s3&
+189._dp*q4*w
C(8,17)=&
-1368._dp*q8&
-2565._dp*q4*p&
+2832._dp*q7*v&
+1728._dp*q3*vp&
-1869._dp*q6*v2&
+384._dp*q5*v3&
+6360._dp*q7*r&
+1728._dp*q3*rp&
-6258._dp*q6*rv&
+1152._dp*q5*rv2&
-4389._dp*q6*r2&
+1152._dp*q5*r2v&
+384._dp*q5*r3&
+2520._dp*q9*s2&
+5184._dp*q5*s2p&
-2268._dp*q8*s2v&
+504._dp*q7*s2v2&
-4536._dp*q8*s2r&
+1008._dp*q7*s2rv&
+504._dp*q7*s2r2&
-16002._dp*q6*sw&
+9396._dp*q5*swv&
+9396._dp*q5*srw&
+2016._dp*q7*s3w&
+2700._dp*q3*w2
C(8,18)=&
+504._dp*q10*s&
+1728._dp*q6*sp&
-504._dp*q9*sv&
+126._dp*q8*sv2&
-1008._dp*q9*sr&
+252._dp*q8*srv&
+126._dp*q8*sr2&
-2286._dp*q7*w&
+1566._dp*q6*wv&
+1566._dp*q6*rw&
+756._dp*q8*s2w
C(8,19)=&
+795._dp*q8&
+432._dp*q4*p&
-894._dp*q7*v&
+192._dp*q6*v2&
-1254._dp*q7*r&
+384._dp*q6*rv&
+192._dp*q6*r2&
-504._dp*q9*s2&
+126._dp*q8*s2v&
+126._dp*q8*s2r&
+1566._dp*q6*sw
C(8,20)=&
-2286._dp*q7*s&
+1566._dp*q6*sv&
+1566._dp*q6*sr&
+252._dp*q8*s3&
+1350._dp*q4*w
C(8,21)=&
+354._dp*q8&
+432._dp*q4*p&
-534._dp*q7*v&
+192._dp*q6*v2&
-894._dp*q7*r&
+384._dp*q6*rv&
+192._dp*q6*r2&
-252._dp*q9*s2&
+126._dp*q8*s2v&
+126._dp*q8*s2r&
+1566._dp*q6*sw
C(8,22)=&
-513._dp*q5&
+432._dp*q4*v&
+432._dp*q4*r&
+864._dp*q6*s2
C(8,23)=&
-88._dp*q10*s&
-2268._dp*q6*sp&
+120._dp*q9*sv&
+1458._dp*q5*svp&
-54._dp*q8*sv2&
+8._dp*q7*sv3&
+1200._dp*q9*sr&
+1458._dp*q5*srp&
-594._dp*q8*srv&
+24._dp*q7*srv2&
-540._dp*q8*sr2&
+24._dp*q7*sr2v&
+8._dp*q7*sr3&
+360._dp*q7*s3p&
+1152._dp*q7*w&
+756._dp*q3*wp&
-1386._dp*q6*wv&
+378._dp*q5*wv2&
-2520._dp*q6*rw&
+756._dp*q5*rwv&
+378._dp*q5*r2w&
-972._dp*q8*s2w&
+432._dp*q7*s2wv&
+432._dp*q7*s2rw&
+1782._dp*q5*sw2
C(8,24)=&
-8._dp*q11&
-324._dp*q7*p&
+12._dp*q10*v&
+243._dp*q6*vp&
-6._dp*q9*v2&
+1._dp*q8*v3&
+120._dp*q10*r&
+243._dp*q6*rp&
-66._dp*q9*rv&
+3._dp*q8*rv2&
-60._dp*q9*r2&
+3._dp*q8*r2v&
+1._dp*q8*r3&
+135._dp*q8*s2p&
-216._dp*q9*sw&
+108._dp*q8*swv&
+108._dp*q8*srw&
+297._dp*q6*w2
C(8,25)=&
+120._dp*q10*s&
+243._dp*q6*sp&
-66._dp*q9*sv&
+3._dp*q8*sv2&
-120._dp*q9*sr&
+6._dp*q8*srv&
+3._dp*q8*sr2&
-360._dp*q7*w&
+126._dp*q6*wv&
+126._dp*q6*rw&
+54._dp*q8*s2w
C(8,26)=&
+144._dp*q8&
+189._dp*q4*p&
-198._dp*q7*v&
+63._dp*q6*v2&
-360._dp*q7*r&
+126._dp*q6*rv&
+63._dp*q6*r2&
-108._dp*q9*s2&
+54._dp*q8*s2v&
+54._dp*q8*s2r&
+594._dp*q6*sw
C(8,27)=&
+12._dp*q10*s&
+243._dp*q6*sp&
-12._dp*q9*sv&
+3._dp*q8*sv2&
-66._dp*q9*sr&
+6._dp*q8*srv&
+3._dp*q8*sr2&
-198._dp*q7*w&
+126._dp*q6*wv&
+126._dp*q6*rw&
+54._dp*q8*s2w
C(8,28)=&
-324._dp*q7*s&
+243._dp*q6*sv&
+243._dp*q6*sr&
+45._dp*q8*s3&
+189._dp*q4*w
C(9,1)=&
+810._dp*q10*sp&
-486._dp*q9*svp&
-486._dp*q9*srp&
-504._dp*q11*w&
-756._dp*q7*wp&
+504._dp*q10*wv&
-126._dp*q9*wv2&
+1098._dp*q10*rw&
-252._dp*q9*rwv&
-126._dp*q9*r2w&
-594._dp*q9*sw2
C(9,2)=&
+16._dp*q12&
+819._dp*q8*p&
+189._dp*q4*p2&
-32._dp*q11*v&
-1008._dp*q7*vp&
+24._dp*q10*v2&
+306._dp*q6*v2p&
-8._dp*q9*v3&
+1._dp*q8*v4&
-500._dp*q11*r&
-1359._dp*q7*rp&
+516._dp*q10*rv&
+612._dp*q6*rvp&
-141._dp*q9*rv2&
+4._dp*q8*rv3&
+654._dp*q10*r2&
+306._dp*q6*r2p&
-258._dp*q9*r2v&
+6._dp*q8*r2v2&
-125._dp*q9*r3&
+4._dp*q8*r3v&
+1._dp*q8*r4&
-1188._dp*q9*s2p&
+675._dp*q8*s2vp&
+675._dp*q8*s2rp&
+972._dp*q10*sw&
+2079._dp*q6*swp&
-972._dp*q9*swv&
+243._dp*q8*swv2&
-1890._dp*q9*srw&
+486._dp*q8*srwv&
+243._dp*q8*sr2w&
-1350._dp*q7*w2&
+675._dp*q6*w2v&
+675._dp*q6*rw2&
+756._dp*q8*s2w2
C(9,3)=&
-1068._dp*q11*s&
-6831._dp*q7*sp&
+2088._dp*q10*sv&
+5265._dp*q6*svp&
-1287._dp*q9*sv2&
+255._dp*q8*sv3&
+5706._dp*q10*sr&
+5265._dp*q6*srp&
-5058._dp*q9*srv&
+765._dp*q8*srv2&
-3771._dp*q9*sr2&
+765._dp*q8*sr2v&
+255._dp*q8*sr3&
+2025._dp*q8*s3p&
+4041._dp*q8*w&
+2835._dp*q4*wp&
-6066._dp*q7*wv&
+2205._dp*q6*wv2&
-8658._dp*q7*rw&
+4410._dp*q6*rwv&
+2205._dp*q6*r2w&
-7290._dp*q9*s2w&
+4320._dp*q8*s2wv&
+4320._dp*q8*s2rw&
+7695._dp*q6*sw2
C(9,4)=&
+8._dp*q12&
+468._dp*q8*p&
+189._dp*q4*p2&
-20._dp*q11*v&
-765._dp*q7*vp&
+18._dp*q10*v2&
+306._dp*q6*v2p&
-7._dp*q9*v3&
+1._dp*q8*v4&
-272._dp*q11*r&
-1116._dp*q7*rp&
+396._dp*q10*rv&
+612._dp*q6*rvp&
-138._dp*q9*rv2&
+4._dp*q8*rv3&
+540._dp*q10*r2&
+306._dp*q6*r2p&
-255._dp*q9*r2v&
+6._dp*q8*r2v2&
-124._dp*q9*r3&
+4._dp*q8*r3v&
+1._dp*q8*r4&
-1053._dp*q9*s2p&
+675._dp*q8*s2vp&
+675._dp*q8*s2rp&
+756._dp*q10*sw&
+2079._dp*q6*swp&
-864._dp*q9*swv&
+243._dp*q8*swv2&
-1782._dp*q9*srw&
+486._dp*q8*srwv&
+243._dp*q8*sr2w&
-1053._dp*q7*w2&
+675._dp*q6*w2v&
+675._dp*q6*rw2&
+756._dp*q8*s2w2
C(9,5)=&
+8100._dp*q9*sp&
-4374._dp*q8*svp&
-4374._dp*q8*srp&
-5544._dp*q10*w&
-5292._dp*q6*wp&
+5040._dp*q9*wv&
-1134._dp*q8*wv2&
+10980._dp*q9*rw&
-2268._dp*q8*rwv&
-1134._dp*q8*r2w&
-5346._dp*q8*sw2
C(9,6)=&
+810._dp*q10*p&
-486._dp*q9*vp&
-486._dp*q9*rp&
-594._dp*q9*w2
C(9,7)=&
-486._dp*q9*sp&
+1098._dp*q10*w&
-252._dp*q9*wv&
-252._dp*q9*rw
C(9,8)=&
-504._dp*q11&
-756._dp*q7*p&
+504._dp*q10*v&
-126._dp*q9*v2&
+1098._dp*q10*r&
-252._dp*q9*rv&
-126._dp*q9*r2&
-1188._dp*q9*sw
C(9,9)=&
-486._dp*q9*sp&
+504._dp*q10*w&
-252._dp*q9*wv&
-252._dp*q9*rw
C(9,10)=&
+810._dp*q10*s&
-486._dp*q9*sv&
-486._dp*q9*sr&
-756._dp*q7*w
C(9,11)=&
+192._dp*q11&
+6552._dp*q7*p&
+756._dp*q3*p2&
-352._dp*q10*v&
-7056._dp*q6*vp&
+240._dp*q9*v2&
+1836._dp*q5*v2p&
-72._dp*q8*v3&
+8._dp*q7*v4&
-5500._dp*q10*r&
-9513._dp*q6*rp&
+5160._dp*q9*rv&
+3672._dp*q5*rvp&
-1269._dp*q8*rv2&
+32._dp*q7*rv3&
+6540._dp*q9*r2&
+1836._dp*q5*r2p&
-2322._dp*q8*r2v&
+48._dp*q7*r2v2&
-1125._dp*q8*r3&
+32._dp*q7*r3v&
+8._dp*q7*r4&
-10692._dp*q8*s2p&
+5400._dp*q7*s2vp&
+5400._dp*q7*s2rp&
+9720._dp*q9*sw&
+12474._dp*q5*swp&
-8748._dp*q8*swv&
+1944._dp*q7*swv2&
-17010._dp*q8*srw&
+3888._dp*q7*srwv&
+1944._dp*q7*sr2w&
-9450._dp*q6*w2&
+4050._dp*q5*w2v&
+4050._dp*q5*rw2&
+6048._dp*q7*s2w2
C(9,12)=&
-2376._dp*q9*sp&
+1350._dp*q8*svp&
+1350._dp*q8*srp&
+972._dp*q10*w&
+2079._dp*q6*wp&
-972._dp*q9*wv&
+243._dp*q8*wv2&
-1890._dp*q9*rw&
+486._dp*q8*rwv&
+243._dp*q8*r2w&
+1512._dp*q8*sw2
C(9,13)=&
-500._dp*q11&
-1359._dp*q7*p&
+516._dp*q10*v&
+612._dp*q6*vp&
-141._dp*q9*v2&
+4._dp*q8*v3&
+1308._dp*q10*r&
+612._dp*q6*rp&
-516._dp*q9*rv&
+12._dp*q8*rv2&
-375._dp*q9*r2&
+12._dp*q8*r2v&
+4._dp*q8*r3&
+675._dp*q8*s2p&
-1890._dp*q9*sw&
+486._dp*q8*swv&
+486._dp*q8*srw&
+675._dp*q6*w2
C(9,14)=&
+972._dp*q10*s&
+2079._dp*q6*sp&
-972._dp*q9*sv&
+243._dp*q8*sv2&
-1890._dp*q9*sr&
+486._dp*q8*srv&
+243._dp*q8*sr2&
-2700._dp*q7*w&
+1350._dp*q6*wv&
+1350._dp*q6*rw&
+1512._dp*q8*s2w
C(9,15)=&
-32._dp*q11&
-1008._dp*q7*p&
+48._dp*q10*v&
+612._dp*q6*vp&
-24._dp*q9*v2&
+4._dp*q8*v3&
+516._dp*q10*r&
+612._dp*q6*rp&
-282._dp*q9*rv&
+12._dp*q8*rv2&
-258._dp*q9*r2&
+12._dp*q8*r2v&
+4._dp*q8*r3&
+675._dp*q8*s2p&
-972._dp*q9*sw&
+486._dp*q8*swv&
+486._dp*q8*srw&
+675._dp*q6*w2
C(9,16)=&
+819._dp*q8&
+378._dp*q4*p&
-1008._dp*q7*v&
+306._dp*q6*v2&
-1359._dp*q7*r&
+612._dp*q6*rv&
+306._dp*q6*r2&
-1188._dp*q9*s2&
+675._dp*q8*s2v&
+675._dp*q8*s2r&
+2079._dp*q6*sw
C(9,17)=&
-11748._dp*q10*s&
-47817._dp*q6*sp&
+20880._dp*q9*sv&
+31590._dp*q5*svp&
-11583._dp*q8*sv2&
+2040._dp*q7*sv3&
+57060._dp*q9*sr&
+31590._dp*q5*srp&
-45522._dp*q8*srv&
+6120._dp*q7*srv2&
-33939._dp*q8*sr2&
+6120._dp*q7*sr2v&
+2040._dp*q7*sr3&
+16200._dp*q7*s3p&
+32328._dp*q7*w&
+11340._dp*q3*wp&
-42462._dp*q6*wv&
+13230._dp*q5*wv2&
-60606._dp*q6*rw&
+26460._dp*q5*rwv&
+13230._dp*q5*r2w&
-65610._dp*q8*s2w&
+34560._dp*q7*s2wv&
+34560._dp*q7*s2rw&
+46170._dp*q5*sw2
C(9,18)=&
-1068._dp*q11&
-6831._dp*q7*p&
+2088._dp*q10*v&
+5265._dp*q6*vp&
-1287._dp*q9*v2&
+255._dp*q8*v3&
+5706._dp*q10*r&
+5265._dp*q6*rp&
-5058._dp*q9*rv&
+765._dp*q8*rv2&
-3771._dp*q9*r2&
+765._dp*q8*r2v&
+255._dp*q8*r3&
+6075._dp*q8*s2p&
-14580._dp*q9*sw&
+8640._dp*q8*swv&
+8640._dp*q8*srw&
+7695._dp*q6*w2
C(9,19)=&
+5706._dp*q10*s&
+5265._dp*q6*sp&
-5058._dp*q9*sv&
+765._dp*q8*sv2&
-7542._dp*q9*sr&
+1530._dp*q8*srv&
+765._dp*q8*sr2&
-8658._dp*q7*w&
+4410._dp*q6*wv&
+4410._dp*q6*rw&
+4320._dp*q8*s2w
C(9,20)=&
+4041._dp*q8&
+2835._dp*q4*p&
-6066._dp*q7*v&
+2205._dp*q6*v2&
-8658._dp*q7*r&
+4410._dp*q6*rv&
+2205._dp*q6*r2&
-7290._dp*q9*s2&
+4320._dp*q8*s2v&
+4320._dp*q8*s2r&
+15390._dp*q6*sw
C(9,21)=&
+2088._dp*q10*s&
+5265._dp*q6*sp&
-2574._dp*q9*sv&
+765._dp*q8*sv2&
-5058._dp*q9*sr&
+1530._dp*q8*srv&
+765._dp*q8*sr2&
-6066._dp*q7*w&
+4410._dp*q6*wv&
+4410._dp*q6*rw&
+4320._dp*q8*s2w
C(9,22)=&
-6831._dp*q7*s&
+5265._dp*q6*sv&
+5265._dp*q6*sr&
+2025._dp*q8*s3&
+2835._dp*q4*w
C(9,23)=&
+96._dp*q11&
+3744._dp*q7*p&
+756._dp*q3*p2&
-220._dp*q10*v&
-5355._dp*q6*vp&
+180._dp*q9*v2&
+1836._dp*q5*v2p&
-63._dp*q8*v3&
+8._dp*q7*v4&
-2992._dp*q10*r&
-7812._dp*q6*rp&
+3960._dp*q9*rv&
+3672._dp*q5*rvp&
-1242._dp*q8*rv2&
+32._dp*q7*rv3&
+5400._dp*q9*r2&
+1836._dp*q5*r2p&
-2295._dp*q8*r2v&
+48._dp*q7*r2v2&
-1116._dp*q8*r3&
+32._dp*q7*r3v&
+8._dp*q7*r4&
-9477._dp*q8*s2p&
+5400._dp*q7*s2vp&
+5400._dp*q7*s2rp&
+7560._dp*q9*sw&
+12474._dp*q5*swp&
-7776._dp*q8*swv&
+1944._dp*q7*swv2&
-16038._dp*q8*srw&
+3888._dp*q7*srwv&
+1944._dp*q7*sr2w&
-7371._dp*q6*w2&
+4050._dp*q5*w2v&
+4050._dp*q5*rw2&
+6048._dp*q7*s2w2
C(9,24)=&
-2106._dp*q9*sp&
+1350._dp*q8*svp&
+1350._dp*q8*srp&
+756._dp*q10*w&
+2079._dp*q6*wp&
-864._dp*q9*wv&
+243._dp*q8*wv2&
-1782._dp*q9*rw&
+486._dp*q8*rwv&
+243._dp*q8*r2w&
+1512._dp*q8*sw2
C(9,25)=&
-272._dp*q11&
-1116._dp*q7*p&
+396._dp*q10*v&
+612._dp*q6*vp&
-138._dp*q9*v2&
+4._dp*q8*v3&
+1080._dp*q10*r&
+612._dp*q6*rp&
-510._dp*q9*rv&
+12._dp*q8*rv2&
-372._dp*q9*r2&
+12._dp*q8*r2v&
+4._dp*q8*r3&
+675._dp*q8*s2p&
-1782._dp*q9*sw&
+486._dp*q8*swv&
+486._dp*q8*srw&
+675._dp*q6*w2
C(9,26)=&
+756._dp*q10*s&
+2079._dp*q6*sp&
-864._dp*q9*sv&
+243._dp*q8*sv2&
-1782._dp*q9*sr&
+486._dp*q8*srv&
+243._dp*q8*sr2&
-2106._dp*q7*w&
+1350._dp*q6*wv&
+1350._dp*q6*rw&
+1512._dp*q8*s2w
C(9,27)=&
-20._dp*q11&
-765._dp*q7*p&
+36._dp*q10*v&
+612._dp*q6*vp&
-21._dp*q9*v2&
+4._dp*q8*v3&
+396._dp*q10*r&
+612._dp*q6*rp&
-276._dp*q9*rv&
+12._dp*q8*rv2&
-255._dp*q9*r2&
+12._dp*q8*r2v&
+4._dp*q8*r3&
+675._dp*q8*s2p&
-864._dp*q9*sw&
+486._dp*q8*swv&
+486._dp*q8*srw&
+675._dp*q6*w2
C(9,28)=&
+468._dp*q8&
+378._dp*q4*p&
-765._dp*q7*v&
+306._dp*q6*v2&
-1116._dp*q7*r&
+612._dp*q6*rv&
+306._dp*q6*r2&
-1053._dp*q9*s2&
+675._dp*q8*s2v&
+675._dp*q8*s2r&
+2079._dp*q6*sw
C(10,1)=&
-16._dp*q15&
-2133._dp*q11*p&
-945._dp*q7*p2&
+32._dp*q14*v&
+2808._dp*q10*vp&
-24._dp*q13*v2&
-918._dp*q9*v2p&
+8._dp*q12*v3&
-1._dp*q11*v4&
+1004._dp*q14*r&
+4509._dp*q10*rp&
-1020._dp*q13*rv&
-1836._dp*q9*rvp&
+267._dp*q12*rv2&
-4._dp*q11*rv3&
-1752._dp*q13*r2&
-918._dp*q9*r2p&
+510._dp*q12*r2v&
-6._dp*q11*r2v2&
+251._dp*q12*r3&
-4._dp*q11*r3v&
-1._dp*q11*r4&
+1350._dp*q12*s2p&
-675._dp*q11*s2vp&
-675._dp*q11*s2rp&
-972._dp*q13*sw&
-6237._dp*q9*swp&
+972._dp*q12*swv&
-243._dp*q11*swv2&
+2484._dp*q12*srw&
-486._dp*q11*srwv&
-243._dp*q11*sr2w&
+4050._dp*q10*w2&
-2025._dp*q9*w2v&
-2025._dp*q9*rw2&
-756._dp*q11*s2w2
C(10,2)=&
+16._dp*q14*s&
+6993._dp*q10*sp&
+3024._dp*q6*sp2&
-32._dp*q13*sv&
-8424._dp*q9*svp&
+24._dp*q12*sv2&
+2511._dp*q8*sv2p&
-8._dp*q11*sv3&
+1._dp*q10*sv4&
-1004._dp*q13*sr&
-12123._dp*q9*srp&
+1020._dp*q12*srv&
+5022._dp*q8*srvp&
-267._dp*q11*srv2&
+4._dp*q10*srv3&
+1752._dp*q12*sr2&
+2511._dp*q8*sr2p&
-510._dp*q11*sr2v&
+6._dp*q10*sr2v2&
-251._dp*q11*sr3&
+4._dp*q10*sr3v&
+1._dp*q10*sr4&
-1350._dp*q11*s3p&
+675._dp*q10*s3vp&
+675._dp*q10*s3rp&
-2040._dp*q11*w&
-8910._dp*q7*wp&
+3060._dp*q10*wv&
+5265._dp*q6*wvp&
-1530._dp*q9*wv2&
+255._dp*q8*wv3&
+9000._dp*q10*rw&
+5265._dp*q6*rwp&
-6030._dp*q9*rwv&
+765._dp*q8*rwv2&
-4500._dp*q9*r2w&
+765._dp*q8*r2wv&
+255._dp*q8*r3w&
+972._dp*q12*s2w&
+9774._dp*q8*s2wp&
-972._dp*q11*s2wv&
+243._dp*q10*s2wv2&
-2484._dp*q11*s2rw&
+486._dp*q10*s2rwv&
+243._dp*q10*s2r2w&
-9990._dp*q9*sw2&
+4995._dp*q8*sw2v&
+4995._dp*q8*srw2&
+756._dp*q10*s3w2&
+2025._dp*q6*w3
C(10,3)=&
+1084._dp*q12&
+11691._dp*q8*p&
+3024._dp*q4*p2&
-3188._dp*q11*v&
-19170._dp*q7*vp&
+3399._dp*q10*v2&
+7776._dp*q6*v2p&
-1550._dp*q9*v3&
+256._dp*q8*v4&
-11315._dp*q11*r&
-24948._dp*q7*rp&
+19434._dp*q10*rv&
+15552._dp*q6*rvp&
-9456._dp*q9*rv2&
+1024._dp*q8*rv3&
+18789._dp*q10*r2&
+7776._dp*q6*r2p&
-14262._dp*q9*r2v&
+1536._dp*q8*r2v2&
-6356._dp*q9*r3&
+1024._dp*q8*r3v&
+256._dp*q8*r4&
-2040._dp*q13*s2&
-38745._dp*q9*s2p&
+3060._dp*q12*s2v&
+26865._dp*q8*s2vp&
-1530._dp*q11*s2v2&
+255._dp*q10*s2v3&
+11700._dp*q12*s2r&
+26865._dp*q8*s2rp&
-7380._dp*q11*s2rv&
+765._dp*q10*s2rv2&
-5850._dp*q11*s2r2&
+765._dp*q10*s2r2v&
+255._dp*q10*s2r3&
+2025._dp*q10*s4p&
+41607._dp*q10*sw&
+47439._dp*q6*swp&
-55350._dp*q9*swv&
+17793._dp*q8*swv2&
-84294._dp*q9*srw&
+35586._dp*q8*srwv&
+17793._dp*q8*sr2w&
-8640._dp*q11*s3w&
+4320._dp*q10*s3wv&
+4320._dp*q10*s3rw&
-29835._dp*q7*w2&
+21600._dp*q6*w2v&
+21600._dp*q6*rw2&
+36801._dp*q8*s2w2
C(10,4)=&
+16._dp*q14*s&
+4860._dp*q10*sp&
+3024._dp*q6*sp2&
-32._dp*q13*sv&
-7074._dp*q9*svp&
+24._dp*q12*sv2&
+2511._dp*q8*sv2p&
-8._dp*q11*sv3&
+1._dp*q10*sv4&
-1004._dp*q13*sr&
-10773._dp*q9*srp&
+1020._dp*q12*srv&
+5022._dp*q8*srvp&
-267._dp*q11*srv2&
+4._dp*q10*srv3&
+1752._dp*q12*sr2&
+2511._dp*q8*sr2p&
-510._dp*q11*sr2v&
+6._dp*q10*sr2v2&
-251._dp*q11*sr3&
+4._dp*q10*sr3v&
+1._dp*q10*sr4&
-1350._dp*q11*s3p&
+675._dp*q10*s3vp&
+675._dp*q10*s3rp&
-1068._dp*q11*w&
-6831._dp*q7*wp&
+2088._dp*q10*wv&
+5265._dp*q6*wvp&
-1287._dp*q9*wv2&
+255._dp*q8*wv3&
+6516._dp*q10*rw&
+5265._dp*q6*rwp&
-5544._dp*q9*rwv&
+765._dp*q8*rwv2&
-4257._dp*q9*r2w&
+765._dp*q8*r2wv&
+255._dp*q8*r3w&
+972._dp*q12*s2w&
+9774._dp*q8*s2wp&
-972._dp*q11*s2wv&
+243._dp*q10*s2wv2&
-2484._dp*q11*s2rw&
+486._dp*q10*s2rwv&
+243._dp*q10*s2r2w&
-8478._dp*q9*sw2&
+4995._dp*q8*sw2v&
+4995._dp*q8*srw2&
+756._dp*q10*s3w2&
+2025._dp*q6*w3
C(10,5)=&
-240._dp*q14&
-23463._dp*q10*p&
-6615._dp*q6*p2&
+448._dp*q13*v&
+28080._dp*q9*vp&
-312._dp*q12*v2&
-8262._dp*q8*v2p&
+96._dp*q11*v3&
-11._dp*q10*v4&
+14056._dp*q13*r&
+45090._dp*q9*rp&
-13260._dp*q12*rv&
-16524._dp*q8*rvp&
+3204._dp*q11*rv2&
-44._dp*q10*rv3&
-22776._dp*q12*r2&
-8262._dp*q8*r2p&
+6120._dp*q11*r2v&
-66._dp*q10*r2v2&
+3012._dp*q11*r3&
-44._dp*q10*r3v&
-11._dp*q10*r4&
+16200._dp*q11*s2p&
-7425._dp*q10*s2vp&
-7425._dp*q10*s2rp&
-12636._dp*q12*sw&
-56133._dp*q8*swp&
+11664._dp*q11*swv&
-2673._dp*q10*swv2&
+29808._dp*q11*srw&
-5346._dp*q10*srwv&
-2673._dp*q10*sr2w&
+40500._dp*q9*w2&
-18225._dp*q8*w2v&
-18225._dp*q8*rw2&
-8316._dp*q10*s2w2
C(10,6)=&
+2700._dp*q12*sp&
-1350._dp*q11*svp&
-1350._dp*q11*srp&
-972._dp*q13*w&
-6237._dp*q9*wp&
+972._dp*q12*wv&
-243._dp*q11*wv2&
+2484._dp*q12*rw&
-486._dp*q11*rwv&
-243._dp*q11*r2w&
-1512._dp*q11*sw2
C(10,7)=&
+1004._dp*q14&
+4509._dp*q10*p&
-1020._dp*q13*v&
-1836._dp*q9*vp&
+267._dp*q12*v2&
-4._dp*q11*v3&
-3504._dp*q13*r&
-1836._dp*q9*rp&
+1020._dp*q12*rv&
-12._dp*q11*rv2&
+753._dp*q12*r2&
-12._dp*q11*r2v&
-4._dp*q11*r3&
-675._dp*q11*s2p&
+2484._dp*q12*sw&
-486._dp*q11*swv&
-486._dp*q11*srw&
-2025._dp*q9*w2
C(10,8)=&
-972._dp*q13*s&
-6237._dp*q9*sp&
+972._dp*q12*sv&
-243._dp*q11*sv2&
+2484._dp*q12*sr&
-486._dp*q11*srv&
-243._dp*q11*sr2&
+8100._dp*q10*w&
-4050._dp*q9*wv&
-4050._dp*q9*rw&
-1512._dp*q11*s2w
C(10,9)=&
+32._dp*q14&
+2808._dp*q10*p&
-48._dp*q13*v&
-1836._dp*q9*vp&
+24._dp*q12*v2&
-4._dp*q11*v3&
-1020._dp*q13*r&
-1836._dp*q9*rp&
+534._dp*q12*rv&
-12._dp*q11*rv2&
+510._dp*q12*r2&
-12._dp*q11*r2v&
-4._dp*q11*r3&
-675._dp*q11*s2p&
+972._dp*q12*sw&
-486._dp*q11*swv&
-486._dp*q11*srw&
-2025._dp*q9*w2
C(10,10)=&
-2133._dp*q11&
-1890._dp*q7*p&
+2808._dp*q10*v&
-918._dp*q9*v2&
+4509._dp*q10*r&
-1836._dp*q9*rv&
-918._dp*q9*r2&
+1350._dp*q12*s2&
-675._dp*q11*s2v&
-675._dp*q11*s2r&
-6237._dp*q9*sw
C(10,11)=&
+224._dp*q13*s&
+69930._dp*q9*sp&
+18144._dp*q5*sp2&
-416._dp*q12*sv&
-75816._dp*q8*svp&
+288._dp*q11*sv2&
+20088._dp*q7*sv2p&
-88._dp*q10*sv3&
+10._dp*q9*sv4&
-13052._dp*q12*sr&
-109107._dp*q8*srp&
+12240._dp*q11*srv&
+40176._dp*q7*srvp&
-2937._dp*q10*srv2&
+40._dp*q9*srv3&
+21024._dp*q11*sr2&
+20088._dp*q7*sr2p&
-5610._dp*q10*sr2v&
+60._dp*q9*sr2v2&
-2761._dp*q10*sr3&
+40._dp*q9*sr3v&
+10._dp*q9*sr4&
-14850._dp*q10*s3p&
+6750._dp*q9*s3vp&
+6750._dp*q9*s3rp&
-22440._dp*q10*w&
-62370._dp*q6*wp&
+30600._dp*q9*wv&
+31590._dp*q5*wvp&
-13770._dp*q8*wv2&
+2040._dp*q7*wv3&
+90000._dp*q9*rw&
+31590._dp*q5*rwp&
-54270._dp*q8*rwv&
+6120._dp*q7*rwv2&
-40500._dp*q8*r2w&
+6120._dp*q7*r2wv&
+2040._dp*q7*r3w&
+11664._dp*q11*s2w&
+78192._dp*q7*s2wp&
-10692._dp*q10*s2wv&
+2430._dp*q9*s2wv2&
-27324._dp*q10*s2rw&
+4860._dp*q9*s2rwv&
+2430._dp*q9*s2r2w&
-89910._dp*q8*sw2&
+39960._dp*q7*sw2v&
+39960._dp*q7*srw2&
+7560._dp*q9*s3w2&
+12150._dp*q5*w3
C(10,12)=&
+16._dp*q14&
+6993._dp*q10*p&
+3024._dp*q6*p2&
-32._dp*q13*v&
-8424._dp*q9*vp&
+24._dp*q12*v2&
+2511._dp*q8*v2p&
-8._dp*q11*v3&
+1._dp*q10*v4&
-1004._dp*q13*r&
-12123._dp*q9*rp&
+1020._dp*q12*rv&
+5022._dp*q8*rvp&
-267._dp*q11*rv2&
+4._dp*q10*rv3&
+1752._dp*q12*r2&
+2511._dp*q8*r2p&
-510._dp*q11*r2v&
+6._dp*q10*r2v2&
-251._dp*q11*r3&
+4._dp*q10*r3v&
+1._dp*q10*r4&
-4050._dp*q11*s2p&
+2025._dp*q10*s2vp&
+2025._dp*q10*s2rp&
+1944._dp*q12*sw&
+19548._dp*q8*swp&
-1944._dp*q11*swv&
+486._dp*q10*swv2&
-4968._dp*q11*srw&
+972._dp*q10*srwv&
+486._dp*q10*sr2w&
-9990._dp*q9*w2&
+4995._dp*q8*w2v&
+4995._dp*q8*rw2&
+2268._dp*q10*s2w2
C(10,13)=&
-1004._dp*q13*s&
-12123._dp*q9*sp&
+1020._dp*q12*sv&
+5022._dp*q8*svp&
-267._dp*q11*sv2&
+4._dp*q10*sv3&
+3504._dp*q12*sr&
+5022._dp*q8*srp&
-1020._dp*q11*srv&
+12._dp*q10*srv2&
-753._dp*q11*sr2&
+12._dp*q10*sr2v&
+4._dp*q10*sr3&
+675._dp*q10*s3p&
+9000._dp*q10*w&
+5265._dp*q6*wp&
-6030._dp*q9*wv&
+765._dp*q8*wv2&
-9000._dp*q9*rw&
+1530._dp*q8*rwv&
+765._dp*q8*r2w&
-2484._dp*q11*s2w&
+486._dp*q10*s2wv&
+486._dp*q10*s2rw&
+4995._dp*q8*sw2
C(10,14)=&
-2040._dp*q11&
-8910._dp*q7*p&
+3060._dp*q10*v&
+5265._dp*q6*vp&
-1530._dp*q9*v2&
+255._dp*q8*v3&
+9000._dp*q10*r&
+5265._dp*q6*rp&
-6030._dp*q9*rv&
+765._dp*q8*rv2&
-4500._dp*q9*r2&
+765._dp*q8*r2v&
+255._dp*q8*r3&
+972._dp*q12*s2&
+9774._dp*q8*s2p&
-972._dp*q11*s2v&
+243._dp*q10*s2v2&
-2484._dp*q11*s2r&
+486._dp*q10*s2rv&
+243._dp*q10*s2r2&
-19980._dp*q9*sw&
+9990._dp*q8*swv&
+9990._dp*q8*srw&
+1512._dp*q10*s3w&
+6075._dp*q6*w2
C(10,15)=&
-32._dp*q13*s&
-8424._dp*q9*sp&
+48._dp*q12*sv&
+5022._dp*q8*svp&
-24._dp*q11*sv2&
+4._dp*q10*sv3&
+1020._dp*q12*sr&
+5022._dp*q8*srp&
-534._dp*q11*srv&
+12._dp*q10*srv2&
-510._dp*q11*sr2&
+12._dp*q10*sr2v&
+4._dp*q10*sr3&
+675._dp*q10*s3p&
+3060._dp*q10*w&
+5265._dp*q6*wp&
-3060._dp*q9*wv&
+765._dp*q8*wv2&
-6030._dp*q9*rw&
+1530._dp*q8*rwv&
+765._dp*q8*r2w&
-972._dp*q11*s2w&
+486._dp*q10*s2wv&
+486._dp*q10*s2rw&
+4995._dp*q8*sw2
C(10,16)=&
+6993._dp*q10*s&
+6048._dp*q6*sp&
-8424._dp*q9*sv&
+2511._dp*q8*sv2&
-12123._dp*q9*sr&
+5022._dp*q8*srv&
+2511._dp*q8*sr2&
-1350._dp*q11*s3&
+675._dp*q10*s3v&
+675._dp*q10*s3r&
-8910._dp*q7*w&
+5265._dp*q6*wv&
+5265._dp*q6*rw&
+9774._dp*q8*s2w
C(10,17)=&
+13008._dp*q11&
+93528._dp*q7*p&
+12096._dp*q3*p2&
-35068._dp*q10*v&
-134190._dp*q6*vp&
+33990._dp*q9*v2&
+46656._dp*q5*v2p&
-13950._dp*q8*v3&
+2048._dp*q7*v4&
-124465._dp*q10*r&
-174636._dp*q6*rp&
+194340._dp*q9*rv&
+93312._dp*q5*rvp&
-85104._dp*q8*rv2&
+8192._dp*q7*rv3&
+187890._dp*q9*r2&
+46656._dp*q5*r2p&
-128358._dp*q8*r2v&
+12288._dp*q7*r2v2&
-57204._dp*q8*r3&
+8192._dp*q7*r3v&
+2048._dp*q7*r4&
-26520._dp*q12*s2&
-348705._dp*q8*s2p&
+36720._dp*q11*s2v&
+214920._dp*q7*s2vp&
-16830._dp*q10*s2v2&
+2550._dp*q9*s2v3&
+140400._dp*q11*s2r&
+214920._dp*q7*s2rp&
-81180._dp*q10*s2rv&
+7650._dp*q9*s2rv2&
-64350._dp*q10*s2r2&
+7650._dp*q9*s2r2v&
+2550._dp*q9*s2r3&
+20250._dp*q9*s4p&
+416070._dp*q9*sw&
+284634._dp*q5*swp&
-498150._dp*q8*swv&
+142344._dp*q7*swv2&
-758646._dp*q8*srw&
+284688._dp*q7*srwv&
+142344._dp*q7*sr2w&
-95040._dp*q10*s3w&
+43200._dp*q9*s3wv&
+43200._dp*q9*s3rw&
-208845._dp*q6*w2&
+129600._dp*q5*w2v&
+129600._dp*q5*rw2&
+294408._dp*q7*s2w2
C(10,18)=&
-4080._dp*q13*s&
-77490._dp*q9*sp&
+6120._dp*q12*sv&
+53730._dp*q8*svp&
-3060._dp*q11*sv2&
+510._dp*q10*sv3&
+23400._dp*q12*sr&
+53730._dp*q8*srp&
-14760._dp*q11*srv&
+1530._dp*q10*srv2&
-11700._dp*q11*sr2&
+1530._dp*q10*sr2v&
+510._dp*q10*sr3&
+8100._dp*q10*s3p&
+41607._dp*q10*w&
+47439._dp*q6*wp&
-55350._dp*q9*wv&
+17793._dp*q8*wv2&
-84294._dp*q9*rw&
+35586._dp*q8*rwv&
+17793._dp*q8*r2w&
-25920._dp*q11*s2w&
+12960._dp*q10*s2wv&
+12960._dp*q10*s2rw&
+73602._dp*q8*sw2
C(10,19)=&
-11315._dp*q11&
-24948._dp*q7*p&
+19434._dp*q10*v&
+15552._dp*q6*vp&
-9456._dp*q9*v2&
+1024._dp*q8*v3&
+37578._dp*q10*r&
+15552._dp*q6*rp&
-28524._dp*q9*rv&
+3072._dp*q8*rv2&
-19068._dp*q9*r2&
+3072._dp*q8*r2v&
+1024._dp*q8*r3&
+11700._dp*q12*s2&
+26865._dp*q8*s2p&
-7380._dp*q11*s2v&
+765._dp*q10*s2v2&
-11700._dp*q11*s2r&
+1530._dp*q10*s2rv&
+765._dp*q10*s2r2&
-84294._dp*q9*sw&
+35586._dp*q8*swv&
+35586._dp*q8*srw&
+4320._dp*q10*s3w&
+21600._dp*q6*w2
C(10,20)=&
+41607._dp*q10*s&
+47439._dp*q6*sp&
-55350._dp*q9*sv&
+17793._dp*q8*sv2&
-84294._dp*q9*sr&
+35586._dp*q8*srv&
+17793._dp*q8*sr2&
-8640._dp*q11*s3&
+4320._dp*q10*s3v&
+4320._dp*q10*s3r&
-59670._dp*q7*w&
+43200._dp*q6*wv&
+43200._dp*q6*rw&
+73602._dp*q8*s2w
C(10,21)=&
-3188._dp*q11&
-19170._dp*q7*p&
+6798._dp*q10*v&
+15552._dp*q6*vp&
-4650._dp*q9*v2&
+1024._dp*q8*v3&
+19434._dp*q10*r&
+15552._dp*q6*rp&
-18912._dp*q9*rv&
+3072._dp*q8*rv2&
-14262._dp*q9*r2&
+3072._dp*q8*r2v&
+1024._dp*q8*r3&
+3060._dp*q12*s2&
+26865._dp*q8*s2p&
-3060._dp*q11*s2v&
+765._dp*q10*s2v2&
-7380._dp*q11*s2r&
+1530._dp*q10*s2rv&
+765._dp*q10*s2r2&
-55350._dp*q9*sw&
+35586._dp*q8*swv&
+35586._dp*q8*srw&
+4320._dp*q10*s3w&
+21600._dp*q6*w2
C(10,22)=&
+11691._dp*q8&
+6048._dp*q4*p&
-19170._dp*q7*v&
+7776._dp*q6*v2&
-24948._dp*q7*r&
+15552._dp*q6*rv&
+7776._dp*q6*r2&
-38745._dp*q9*s2&
+26865._dp*q8*s2v&
+26865._dp*q8*s2r&
+2025._dp*q10*s4&
+47439._dp*q6*sw
C(10,23)=&
+224._dp*q13*s&
+48600._dp*q9*sp&
+18144._dp*q5*sp2&
-416._dp*q12*sv&
-63666._dp*q8*svp&
+288._dp*q11*sv2&
+20088._dp*q7*sv2p&
-88._dp*q10*sv3&
+10._dp*q9*sv4&
-13052._dp*q12*sr&
-96957._dp*q8*srp&
+12240._dp*q11*srv&
+40176._dp*q7*srvp&
-2937._dp*q10*srv2&
+40._dp*q9*srv3&
+21024._dp*q11*sr2&
+20088._dp*q7*sr2p&
-5610._dp*q10*sr2v&
+60._dp*q9*sr2v2&
-2761._dp*q10*sr3&
+40._dp*q9*sr3v&
+10._dp*q9*sr4&
-14850._dp*q10*s3p&
+6750._dp*q9*s3vp&
+6750._dp*q9*s3rp&
-11748._dp*q10*w&
-47817._dp*q6*wp&
+20880._dp*q9*wv&
+31590._dp*q5*wvp&
-11583._dp*q8*wv2&
+2040._dp*q7*wv3&
+65160._dp*q9*rw&
+31590._dp*q5*rwp&
-49896._dp*q8*rwv&
+6120._dp*q7*rwv2&
-38313._dp*q8*r2w&
+6120._dp*q7*r2wv&
+2040._dp*q7*r3w&
+11664._dp*q11*s2w&
+78192._dp*q7*s2wp&
-10692._dp*q10*s2wv&
+2430._dp*q9*s2wv2&
-27324._dp*q10*s2rw&
+4860._dp*q9*s2rwv&
+2430._dp*q9*s2r2w&
-76302._dp*q8*sw2&
+39960._dp*q7*sw2v&
+39960._dp*q7*srw2&
+7560._dp*q9*s3w2&
+12150._dp*q5*w3
C(10,24)=&
+16._dp*q14&
+4860._dp*q10*p&
+3024._dp*q6*p2&
-32._dp*q13*v&
-7074._dp*q9*vp&
+24._dp*q12*v2&
+2511._dp*q8*v2p&
-8._dp*q11*v3&
+1._dp*q10*v4&
-1004._dp*q13*r&
-10773._dp*q9*rp&
+1020._dp*q12*rv&
+5022._dp*q8*rvp&
-267._dp*q11*rv2&
+4._dp*q10*rv3&
+1752._dp*q12*r2&
+2511._dp*q8*r2p&
-510._dp*q11*r2v&
+6._dp*q10*r2v2&
-251._dp*q11*r3&
+4._dp*q10*r3v&
+1._dp*q10*r4&
-4050._dp*q11*s2p&
+2025._dp*q10*s2vp&
+2025._dp*q10*s2rp&
+1944._dp*q12*sw&
+19548._dp*q8*swp&
-1944._dp*q11*swv&
+486._dp*q10*swv2&
-4968._dp*q11*srw&
+972._dp*q10*srwv&
+486._dp*q10*sr2w&
-8478._dp*q9*w2&
+4995._dp*q8*w2v&
+4995._dp*q8*rw2&
+2268._dp*q10*s2w2
C(10,25)=&
-1004._dp*q13*s&
-10773._dp*q9*sp&
+1020._dp*q12*sv&
+5022._dp*q8*svp&
-267._dp*q11*sv2&
+4._dp*q10*sv3&
+3504._dp*q12*sr&
+5022._dp*q8*srp&
-1020._dp*q11*srv&
+12._dp*q10*srv2&
-753._dp*q11*sr2&
+12._dp*q10*sr2v&
+4._dp*q10*sr3&
+675._dp*q10*s3p&
+6516._dp*q10*w&
+5265._dp*q6*wp&
-5544._dp*q9*wv&
+765._dp*q8*wv2&
-8514._dp*q9*rw&
+1530._dp*q8*rwv&
+765._dp*q8*r2w&
-2484._dp*q11*s2w&
+486._dp*q10*s2wv&
+486._dp*q10*s2rw&
+4995._dp*q8*sw2
C(10,26)=&
-1068._dp*q11&
-6831._dp*q7*p&
+2088._dp*q10*v&
+5265._dp*q6*vp&
-1287._dp*q9*v2&
+255._dp*q8*v3&
+6516._dp*q10*r&
+5265._dp*q6*rp&
-5544._dp*q9*rv&
+765._dp*q8*rv2&
-4257._dp*q9*r2&
+765._dp*q8*r2v&
+255._dp*q8*r3&
+972._dp*q12*s2&
+9774._dp*q8*s2p&
-972._dp*q11*s2v&
+243._dp*q10*s2v2&
-2484._dp*q11*s2r&
+486._dp*q10*s2rv&
+243._dp*q10*s2r2&
-16956._dp*q9*sw&
+9990._dp*q8*swv&
+9990._dp*q8*srw&
+1512._dp*q10*s3w&
+6075._dp*q6*w2
C(10,27)=&
-32._dp*q13*s&
-7074._dp*q9*sp&
+48._dp*q12*sv&
+5022._dp*q8*svp&
-24._dp*q11*sv2&
+4._dp*q10*sv3&
+1020._dp*q12*sr&
+5022._dp*q8*srp&
-534._dp*q11*srv&
+12._dp*q10*srv2&
-510._dp*q11*sr2&
+12._dp*q10*sr2v&
+4._dp*q10*sr3&
+675._dp*q10*s3p&
+2088._dp*q10*w&
+5265._dp*q6*wp&
-2574._dp*q9*wv&
+765._dp*q8*wv2&
-5544._dp*q9*rw&
+1530._dp*q8*rwv&
+765._dp*q8*r2w&
-972._dp*q11*s2w&
+486._dp*q10*s2wv&
+486._dp*q10*s2rw&
+4995._dp*q8*sw2
C(10,28)=&
+4860._dp*q10*s&
+6048._dp*q6*sp&
-7074._dp*q9*sv&
+2511._dp*q8*sv2&
-10773._dp*q9*sr&
+5022._dp*q8*srv&
+2511._dp*q8*sr2&
-1350._dp*q11*s3&
+675._dp*q10*s3v&
+675._dp*q10*s3r&
-6831._dp*q7*w&
+5265._dp*q6*wv&
+5265._dp*q6*rw&
+9774._dp*q8*s2w
end subroutine fg10ex