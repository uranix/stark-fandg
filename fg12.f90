subroutine fg12(q,s,r,w,v,p,C)
implicit double precision (q,s,r,w,v,p)
integer, parameter :: dp = kind(1.d0)
real(dp), intent(in) :: q,s,r,w,v,p
real(dp), intent(out) :: C(0:12,4)
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
q16=q15*q
q17=q16*q
q18=q17*q
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
vp2=vp*p
r3=r2*r
r2w=r2*w
r2p=r2*p
rw2=rw*w
rwp=rw*p
rp2=rp*p
w3=w2*w
w2p=w2*p
wp2=wp*p
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
sw2v=swv*w
swvp=swv*p
svp2=svp*p
sr3=sr2*r
s2r2=sr2*s
sr2w=sr2*w
sr2p=sr2*p
s2rp=s2r*p
srw2=srw*w
srwp=srw*p
srp2=srp*p
s2wv=s2w*v
s2rw=s2w*r
s3w=s2w*s
s2w2=s2w*w
s2wp=s2w*p
s3p=s2p*s
s2p2=s2p*p
sw3=sw2*w
sw2p=sw2*p
v4=v3*v
rv3=v3*r
wv3=v3*w
v3p=v3*p
r2v2=rv2*r
rwv2=rv2*w
rv2p=rv2*p
w2v2=wv2*w
wv2p=wv2*p
r3v=r2v*r
r2wv=r2v*w
r2vp=r2v*p
rw2v=rwv*w
rwvp=rwv*p
w3v=w2v*w
r4=r3*r
r3w=r3*w
r3p=r3*p
r2w2=r2w*w
r2wp=r2w*p
rw3=rw2*w
sv4=sv3*v
srv3=sv3*r
s2v3=sv3*s
swv3=sv3*w
sv3p=sv3*p
sr2v2=srv2*r
s2rv2=srv2*s
srwv2=srv2*w
srv2p=srv2*p
s2v2p=s2v2*p
sw2v2=swv2*w
sr3v=sr2v*r
s2r2v=sr2v*s
sr2wv=sr2v*w
sr2vp=sr2v*p
s2rvp=s2rv*p
srw2v=srwv*w
sr4=sr3*r
s2r3=sr3*s
sr3w=sr3*w
sr3p=sr3*p
s2r2p=s2r2*p
sr2w2=sr2w*w
s2wv2=s2wv*v
s2rwv=s2wv*r
s3wv=s2wv*s
s2w2v=s2wv*w
s2wvp=s2wv*p
s2r2w=s2rw*r
s3rw=s2rw*s
s2rw2=s2rw*w
s2rwp=s2rw*p
s3w2=s3w*w
s2w3=s2w2*w
s3vp=s3p*v
s3rp=s3p*r
s4p=s3p*s
s3wp=s3p*w
s3p2=s3p*p
v5=v4*v
rv4=v4*r
wv4=v4*w
r2v3=rv3*r
rwv3=rv3*w
r3v2=r2v2*r
r2wv2=r2v2*w
r4v=r3v*r
r3wv=r3v*w
r5=r4*r
r4w=r4*w
sv5=sv4*v
srv4=sv4*r
s2v4=sv4*s
sr2v3=srv3*r
s2rv3=srv3*s
sr3v2=sr2v2*r
s2r2v2=sr2v2*s
sr4v=sr3v*r
s2r3v=sr3v*s
sr5=sr4*r
s2r4=sr4*s
s2wv3=s2wv2*v
s2rwv2=s2wv2*r
s3wv2=s2wv2*s
s2r2wv=s2rwv*r
s3rwv=s2rwv*s
s2r3w=s2r2w*r
s3r2w=s2r2w*s
s3w2v=s3w2*v
s3rw2=s3w2*r
s4w2=s3w2*s
s3v2p=s3vp*v
s3rvp=s3vp*r
s4vp=s3vp*s
s3r2p=s3rp*r
s4rp=s3rp*s
s4wp=s4p*w
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
C(11,1)=&
-14148._dp*q13*sp&
-12096._dp*q9*sp2&
+17118._dp*q12*svp&
-5022._dp*q11*sv2p&
+28566._dp*q12*srp&
-10044._dp*q11*srvp&
-5022._dp*q11*sr2p&
+4080._dp*q14*w&
+33372._dp*q10*wp&
-6120._dp*q13*wv&
-21060._dp*q9*wvp&
+3060._dp*q12*wv2&
-510._dp*q11*wv3&
-26100._dp*q13*rw&
-21060._dp*q9*rwp&
+16110._dp*q12*rwv&
-1530._dp*q11*rwv2&
+13050._dp*q12*r2w&
-1530._dp*q11*r2wv&
-510._dp*q11*r3w&
-19548._dp*q11*s2wp&
+19980._dp*q12*sw2&
-9990._dp*q11*sw2v&
-9990._dp*q11*srw2&
-8100._dp*q9*w3
C(11,2)=&
-32._dp*q15&
-11166._dp*q11*p&
-12879._dp*q7*p2&
+80._dp*q14*v&
+21285._dp*q10*vp&
+8289._dp*q6*vp2&
-80._dp*q13*v2&
-13383._dp*q9*v2p&
+40._dp*q12*v3&
+2766._dp*q8*v3p&
-10._dp*q11*v4&
+1._dp*q10*v5&
+4064._dp*q14*r&
+41535._dp*q10*rp&
+8289._dp*q6*rp2&
-6136._dp*q13*rv&
-38700._dp*q9*rvp&
+3108._dp*q12*rv2&
+8298._dp*q8*rv2p&
-538._dp*q11*rv3&
+5._dp*q10*rv4&
-13508._dp*q13*r2&
-25317._dp*q9*r2p&
+9822._dp*q12*r2v&
+8298._dp*q8*r2vp&
-1554._dp*q11*r2v2&
+10._dp*q10*r2v3&
+6754._dp*q12*r3&
+2766._dp*q8*r3p&
-1534._dp*q11*r3v&
+10._dp*q10*r3v2&
-508._dp*q11*r4&
+5._dp*q10*r4v&
+1._dp*q10*r5&
+23220._dp*q12*s2p&
+21870._dp*q8*s2p2&
-26190._dp*q11*s2vp&
+7290._dp*q10*s2v2p&
-41175._dp*q11*s2rp&
+14580._dp*q10*s2rvp&
+7290._dp*q10*s2r2p&
-8064._dp*q13*sw&
-111996._dp*q9*swp&
+12096._dp*q12*swv&
+65664._dp*q8*swvp&
-6048._dp*q11*swv2&
+1008._dp*q10*swv3&
+46980._dp*q12*srw&
+65664._dp*q8*srwp&
-29538._dp*q11*srwv&
+3024._dp*q10*srwv2&
-23490._dp*q11*sr2w&
+3024._dp*q10*sr2wv&
+1008._dp*q10*sr3w&
+23085._dp*q10*s3wp&
+29160._dp*q10*w2&
+21870._dp*q6*w2p&
-29160._dp*q9*w2v&
+7290._dp*q8*w2v2&
-44145._dp*q9*rw2&
+14580._dp*q8*rw2v&
+7290._dp*q8*r2w2&
-27432._dp*q11*s2w2&
+13716._dp*q10*s2w2v&
+13716._dp*q10*s2rw2&
+23085._dp*q8*sw3
C(11,3)=&
+8304._dp*q14*s&
+187812._dp*q10*sp&
+68607._dp*q6*sp2&
-20640._dp*q13*sv&
-283608._dp*q9*svp&
+18504._dp*q12*sv2&
+105138._dp*q8*sv2p&
-7176._dp*q11*sv3&
+1023._dp*q10*sv4&
-97644._dp*q13*sr&
-381024._dp*q9*srp&
+138474._dp*q12*srv&
+210276._dp*q8*srvp&
-53010._dp*q11*srv2&
+4092._dp*q10*srv3&
+149670._dp*q12*sr2&
+105138._dp*q8*sr2p&
-84492._dp*q11*sr2v&
+6138._dp*q10*sr2v2&
-38658._dp*q11*sr3&
+4092._dp*q10*sr3v&
+1023._dp*q10*sr4&
-110565._dp*q11*s3p&
+66825._dp*q10*s3vp&
+66825._dp*q10*s3rp&
-61338._dp*q11*w&
-179307._dp*q7*wp&
+133047._dp*q10*wv&
+142560._dp*q6*wvp&
-93429._dp*q9*wv2&
+21120._dp*q8*wv3&
+271017._dp*q10*rw&
+142560._dp*q6*rwp&
-290808._dp*q9*rwv&
+63360._dp*q8*rwv2&
-197379._dp*q9*r2w&
+63360._dp*q8*r2wv&
+21120._dp*q8*r3w&
+155412._dp*q12*s2w&
+353727._dp*q8*s2wp&
-179874._dp*q11*s2wv&
+51084._dp*q10*s2wv2&
-296892._dp*q11*s2rw&
+102168._dp*q10*s2rwv&
+51084._dp*q10*s2r2w&
-411291._dp*q9*sw2&
+271755._dp*q8*sw2v&
+271755._dp*q8*srw2&
+87318._dp*q10*s3w2&
+66825._dp*q6*w3
C(11,4)=&
-16._dp*q15&
-5928._dp*q11*p&
-9855._dp*q7*p2&
+48._dp*q14*v&
+14022._dp*q10*vp&
+8289._dp*q6*vp2&
-56._dp*q13*v2&
-10872._dp*q9*v2p&
+32._dp*q12*v3&
+2766._dp*q8*v3p&
-9._dp*q11*v4&
+1._dp*q10*v5&
+2088._dp*q14*r&
+28980._dp*q10*rp&
+8289._dp*q6*rp2&
-4144._dp*q13*rv&
-33678._dp*q9*rvp&
+2598._dp*q12*rv2&
+8298._dp*q8*rv2p&
-534._dp*q11*rv3&
+5._dp*q10*rv4&
-9272._dp*q13*r2&
-22806._dp*q9*r2p&
+8826._dp*q12*r2v&
+8298._dp*q8*r2vp&
-1548._dp*q11*r2v2&
+10._dp*q10*r2v3&
+6260._dp*q12*r3&
+2766._dp*q8*r3p&
-1530._dp*q11*r3v&
+10._dp*q10*r3v2&
-507._dp*q11*r4&
+5._dp*q10*r4v&
+1._dp*q10*r5&
+19170._dp*q12*s2p&
+21870._dp*q8*s2p2&
-24165._dp*q11*s2vp&
+7290._dp*q10*s2v2p&
-39150._dp*q11*s2rp&
+14580._dp*q10*s2rvp&
+7290._dp*q10*s2r2p&
-6120._dp*q13*sw&
-92448._dp*q9*swp&
+10152._dp*q12*swv&
+65664._dp*q8*swvp&
-5562._dp*q11*swv2&
+1008._dp*q10*swv3&
+40500._dp*q12*srw&
+65664._dp*q8*srwp&
-28566._dp*q11*srwv&
+3024._dp*q10*srwv2&
-23004._dp*q11*sr2w&
+3024._dp*q10*sr2wv&
+1008._dp*q10*sr3w&
+23085._dp*q10*s3wp&
+19170._dp*q10*w2&
+21870._dp*q6*w2p&
-24165._dp*q9*w2v&
+7290._dp*q8*w2v2&
-39150._dp*q9*rw2&
+14580._dp*q8*rw2v&
+7290._dp*q8*r2w2&
-25164._dp*q11*s2w2&
+13716._dp*q10*s2w2v&
+13716._dp*q10*s2rw2&
+23085._dp*q8*sw3
C(12,1)=&
+32._dp*q18&
+29394._dp*q14*p&
+58347._dp*q10*p2&
-80._dp*q17*v&
-58671._dp*q13*vp&
-41445._dp*q9*vp2&
+80._dp*q16*v2&
+38583._dp*q12*v2p&
-40._dp*q15*v3&
-8298._dp*q11*v3p&
+10._dp*q14*v4&
-1._dp*q13*v5&
-8144._dp*q17*r&
-143721._dp*q13*rp&
-41445._dp*q9*rp2&
+12256._dp*q16*rv&
+131598._dp*q12*rvp&
-6168._dp*q15*rv2&
-24894._dp*q11*rv2p&
+1048._dp*q14*rv3&
-5._dp*q13*rv4&
+39608._dp*q16*r2&
+93015._dp*q12*r2p&
-25932._dp*q15*r2v&
-24894._dp*q11*r2vp&
+3084._dp*q14*r2v2&
-10._dp*q13*r2v3&
-19804._dp*q15*r3&
-8298._dp*q11*r3p&
+3064._dp*q14*r3v&
-10._dp*q13*r3v2&
+1018._dp*q14*r4&
-5._dp*q13*r4v&
-1._dp*q13*r5&
-29160._dp*q15*s2p&
-65610._dp*q11*s2p2&
+29160._dp*q14*s2vp&
-7290._dp*q13*s2v2p&
+52245._dp*q14*s2rp&
-14580._dp*q13*s2rvp&
-7290._dp*q13*s2r2p&
+8064._dp*q16*sw&
+329346._dp*q12*swp&
-12096._dp*q15*swv&
-196992._dp*q11*swvp&
+6048._dp*q14*swv2&
-1008._dp*q13*swv3&
-66960._dp*q15*srw&
-196992._dp*q11*srwp&
+39528._dp*q14*srwv&
-3024._dp*q13*srwv2&
+33480._dp*q14*sr2w&
-3024._dp*q13*sr2wv&
-1008._dp*q13*sr3w&
-23085._dp*q13*s3wp&
-87480._dp*q13*w2&
-109350._dp*q9*w2p&
+87480._dp*q12*w2v&
-21870._dp*q11*w2v2&
+156735._dp*q12*rw2&
-43740._dp*q11*rw2v&
-21870._dp*q11*r2w2&
+27432._dp*q14*s2w2&
-13716._dp*q13*s2w2v&
-13716._dp*q13*s2rw2&
-69255._dp*q11*sw3
C(12,2)=&
-32._dp*q17*s&
-111222._dp*q13*sp&
-235926._dp*q9*sp2&
+80._dp*q16*sv&
+202851._dp*q12*svp&
+150849._dp*q8*svp2&
-80._dp*q15*sv2&
-121392._dp*q11*sv2p&
+40._dp*q14*sv3&
+23886._dp*q10*sv3p&
-10._dp*q13*sv4&
+1._dp*q12*sv5&
+8144._dp*q16*sr&
+435267._dp*q12*srp&
+150849._dp*q8*srp2&
-12256._dp*q15*srv&
-379728._dp*q11*srvp&
+6168._dp*q14*srv2&
+71658._dp*q10*srv2p&
-1048._dp*q13*srv3&
+5._dp*q12*srv4&
-39608._dp*q15*sr2&
-258336._dp*q11*sr2p&
+25932._dp*q14*sr2v&
+71658._dp*q10*sr2vp&
-3084._dp*q13*sr2v2&
+10._dp*q12*sr2v3&
+19804._dp*q14*sr3&
+23886._dp*q10*sr3p&
-3064._dp*q13*sr3v&
+10._dp*q12*sr3v2&
-1018._dp*q13*sr4&
+5._dp*q12*sr4v&
+1._dp*q12*sr5&
+29160._dp*q14*s3p&
+88695._dp*q10*s3p2&
-29160._dp*q13*s3vp&
+7290._dp*q12*s3v2p&
-52245._dp*q13*s3rp&
+14580._dp*q12*s3rvp&
+7290._dp*q12*s3r2p&
+16368._dp*q14*w&
+287793._dp*q10*wp&
+68607._dp*q6*wp2&
-32736._dp*q13*wv&
-349272._dp*q9*wvp&
+24552._dp*q12*wv2&
+105138._dp*q8*wv2p&
-8184._dp*q11*wv3&
+1023._dp*q10*wv4&
-178752._dp*q13*rw&
-458784._dp*q9*rwp&
+195120._dp*q12*rwv&
+210276._dp*q8*rwvp&
-61056._dp*q11*rwv2&
+4092._dp*q10*rwv3&
+211716._dp*q12*r2w&
+105138._dp*q8*r2wp&
-97560._dp*q11*r2wv&
+6138._dp*q10*r2wv2&
-44688._dp*q11*r3w&
+4092._dp*q10*r3wv&
+1023._dp*q10*r4w&
-8064._dp*q15*s2w&
-592542._dp*q11*s2wp&
+12096._dp*q14*s2wv&
+337419._dp*q10*s2wvp&
-6048._dp*q13*s2wv2&
+1008._dp*q12*s2wv3&
+66960._dp*q14*s2rw&
+337419._dp*q10*s2rwp&
-39528._dp*q13*s2rwv&
+3024._dp*q12*s2rwv2&
-33480._dp*q13*s2r2w&
+3024._dp*q12*s2r2wv&
+1008._dp*q12*s2r3w&
+23085._dp*q12*s4wp&
+233496._dp*q12*sw2&
+375597._dp*q8*sw2p&
-233496._dp*q11*sw2v&
+58374._dp*q10*sw2v2&
-385047._dp*q11*srw2&
+116748._dp*q10*srw2v&
+58374._dp*q10*sr2w2&
-27432._dp*q13*s3w2&
+13716._dp*q12*s3w2v&
+13716._dp*q12*s3rw2&
-133650._dp*q9*w3&
+66825._dp*q8*w3v&
+66825._dp*q8*rw3&
+110403._dp*q10*s2w3
C(12,3)=&
-8336._dp*q15&
-260316._dp*q11*p&
-260793._dp*q7*p2&
+29024._dp*q14*v&
+625752._dp*q10*vp&
+219456._dp*q6*vp2&
-39224._dp*q13*v2&
-495558._dp*q9*v2p&
+25720._dp*q12*v3&
+129024._dp*q8*v3p&
-8209._dp*q11*v4&
+1024._dp*q10*v5&
+171350._dp*q14*r&
+1.060695e6_dp*q10*rp&
+219456._dp*q6*rp2&
-395941._dp*q13*rv&
-1.346976e6_dp*q9*rvp&
+306525._dp*q12*rv2&
+387072._dp*q8*rv2p&
-85936._dp*q11*rv3&
+5120._dp*q10*rv4&
-531839._dp*q13*r2&
-851418._dp*q9*r2p&
+673266._dp*q12*r2v&
+387072._dp*q8*r2vp&
-208554._dp*q11*r2v2&
+10240._dp*q10*r2v3&
+392461._dp*q12*r3&
+129024._dp*q8*r3p&
-192136._dp*q11*r3v&
+10240._dp*q10*r3v2&
-61309._dp*q11*r4&
+5120._dp*q10*r4v&
+1024._dp*q10*r5&
+16368._dp*q16*s2&
+1.265355e6_dp*q12*s2p&
+718632._dp*q8*s2p2&
-32736._dp*q15*s2v&
-1.726002e6_dp*q11*s2vp&
+24552._dp*q14*s2v2&
+574263._dp*q10*s2v2p&
-8184._dp*q13*s2v3&
+1023._dp*q12*s2v4&
-237072._dp*q15*s2r&
-2.406564e6_dp*q11*s2rp&
+253440._dp*q14*s2rv&
+1.148526e6_dp*q10*s2rvp&
-75636._dp*q13*s2rv2&
+4092._dp*q12*s2rv3&
+316206._dp*q14*s2r2&
+574263._dp*q10*s2r2p&
-126720._dp*q13*s2r2v&
+6138._dp*q12*s2r2v2&
-59268._dp*q13*s2r3&
+4092._dp*q12*s2r3v&
+1023._dp*q12*s2r4&
-133650._dp*q13*s4p&
+66825._dp*q12*s4vp&
+66825._dp*q12*s4rp&
-723906._dp*q13*sw&
-3.59262e6_dp*q9*swp&
+1.401921e6_dp*q12*swv&
+2.660256e6_dp*q8*swvp&
-877608._dp*q11*swv2&
+178812._dp*q10*swv3&
+3.203091e6_dp*q12*srw&
+2.660256e6_dp*q8*srwp&
-2.946996e6_dp*q11*srwv&
+536436._dp*q10*srwv2&
-2.069388e6_dp*q11*sr2w&
+536436._dp*q10*sr2wv&
+178812._dp*q10*sr3w&
+204336._dp*q14*s3w&
+1.459377e6_dp*q10*s3wp&
-204336._dp*q13*s3wv&
+51084._dp*q12*s3wv2&
-378972._dp*q13*s3rw&
+102168._dp*q12*s3rwv&
+51084._dp*q12*s3r2w&
+977562._dp*q10*w2&
+650025._dp*q6*w2p&
-1.37673e6_dp*q9*w2v&
+469125._dp*q8*w2v2&
-1.90404e6_dp*q9*rw2&
+938250._dp*q8*rw2v&
+469125._dp*q8*r2w2&
-2.312118e6_dp*q11*s2w2&
+1.397439e6_dp*q10*s2w2v&
+1.397439e6_dp*q10*s2rw2&
+87318._dp*q12*s4w2&
+1.172475e6_dp*q8*sw3
C(12,4)=&
-32._dp*q17*s&
-72504._dp*q13*sp&
-192186._dp*q9*sp2&
+80._dp*q16*sv&
+154332._dp*q12*svp&
+150849._dp*q8*svp2&
-80._dp*q15*sv2&
-106812._dp*q11*sv2p&
+40._dp*q14*sv3&
+23886._dp*q10*sv3p&
-10._dp*q13*sv4&
+1._dp*q12*sv5&
+8144._dp*q16*sr&
+345924._dp*q12*srp&
+150849._dp*q8*srp2&
-12256._dp*q15*srv&
-350568._dp*q11*srvp&
+6168._dp*q14*srv2&
+71658._dp*q10*srv2p&
-1048._dp*q13*srv3&
+5._dp*q12*srv4&
-39608._dp*q15*sr2&
-243756._dp*q11*sr2p&
+25932._dp*q14*sr2v&
+71658._dp*q10*sr2vp&
-3084._dp*q13*sr2v2&
+10._dp*q12*sr2v3&
+19804._dp*q14*sr3&
+23886._dp*q10*sr3p&
-3064._dp*q13*sr3v&
+10._dp*q12*sr3v2&
-1018._dp*q13*sr4&
+5._dp*q12*sr4v&
+1._dp*q12*sr5&
+29160._dp*q14*s3p&
+88695._dp*q10*s3p2&
-29160._dp*q13*s3vp&
+7290._dp*q12*s3v2p&
-52245._dp*q13*s3rp&
+14580._dp*q12*s3rvp&
+7290._dp*q12*s3r2p&
+8304._dp*q14*w&
+187812._dp*q10*wp&
+68607._dp*q6*wp2&
-20640._dp*q13*wv&
-283608._dp*q9*wvp&
+18504._dp*q12*wv2&
+105138._dp*q8*wv2p&
-7176._dp*q11*wv3&
+1023._dp*q10*wv4&
-111792._dp*q13*rw&
-393120._dp*q9*rwp&
+155592._dp*q12*rwv&
+210276._dp*q8*rwvp&
-58032._dp*q11*rwv2&
+4092._dp*q10*rwv3&
+178236._dp*q12*r2w&
+105138._dp*q8*r2wp&
-94536._dp*q11*r2wv&
+6138._dp*q10*r2wv2&
-43680._dp*q11*r3w&
+4092._dp*q10*r3wv&
+1023._dp*q10*r4w&
-8064._dp*q15*s2w&
-523287._dp*q11*s2wp&
+12096._dp*q14*s2wv&
+337419._dp*q10*s2wvp&
-6048._dp*q13*s2wv2&
+1008._dp*q12*s2wv3&
+66960._dp*q14*s2rw&
+337419._dp*q10*s2rwp&
-39528._dp*q13*s2rwv&
+3024._dp*q12*s2rwv2&
-33480._dp*q13*s2r2w&
+3024._dp*q12*s2r2wv&
+1008._dp*q12*s2r3w&
+23085._dp*q12*s4wp&
+178632._dp*q12*sw2&
+375597._dp*q8*sw2p&
-206064._dp*q11*sw2v&
+58374._dp*q10*sw2v2&
-357615._dp*q11*srw2&
+116748._dp*q10*srw2v&
+58374._dp*q10*sr2w2&
-27432._dp*q13*s3w2&
+13716._dp*q12*s3w2v&
+13716._dp*q12*s3rw2&
-110565._dp*q9*w3&
+66825._dp*q8*w3v&
+66825._dp*q8*rw3&
+110403._dp*q10*s2w3
end subroutine fg12
