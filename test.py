from __future__ import print_function
from stark import fandg, fandgex
import numpy as np

def main():
    r0 = [1, 2, 3]
    v0 = [-1, 1, 2]
    p = [1, 1, 1]
    tau = 0.01
    alpha = 2
    kmax = 18
    delta = 1e-6
    dvec = np.array([-0.470387, -0.9156114, 0.71380349])

    r, v, t, drdtau, dvdtau, dtdtau, \
        drdr0, dvdr0, dtdr0, \
        drdv0, dvdv0, dtdv0, \
        drdp, dvdp, dtdp, im \
        = fandgex(r0, v0, p, tau, kmax, alpha)

    print('r =', r)
    print('v =', v)
    print('t =', t)
    print('|H0 - H| / |H0| =', abs(im[0,0] - im[0,1]) / abs(im[0,0]))

    r1,v1,t1,_ = fandg(r0, v0, p, tau * (1 - delta), kmax, alpha)
    r2,v2,t2,_ = fandg(r0, v0, p, tau * (1 + delta), kmax, alpha)

    print('drdtau ~', ((r2-r1)/(2*delta) - tau*drdtau) / drdtau)
    print('dvdtau ~', ((v2-v1)/(2*delta) - tau*dvdtau) / dvdtau)
    print('dtdtau  ~', ((t2-t1)/(2*delta) - tau*dtdtau) / dtdtau)

    r1,v1,t1,_ = fandg(r0 - dvec * delta, v0, p, tau, kmax, alpha)
    r2,v2,t2,_ = fandg(r0 + dvec * delta, v0, p, tau, kmax, alpha)

    print('drdr0 ~', ((r2-r1)/(2*delta) - drdr0.dot(dvec)) / drdr0.dot(dvec))
    print('dvdr0 ~', ((v2-v1)/(2*delta) - dvdr0.dot(dvec)) / dvdr0.dot(dvec))
    print('dtdr0 ~', ((t2-t1)/(2*delta) - dtdr0.dot(dvec)) / dtdr0.dot(dvec))

    r1,v1,t1,_ = fandg(r0, v0 - dvec * delta, p, tau, kmax, alpha)
    r2,v2,t2,_ = fandg(r0, v0 + dvec * delta, p, tau, kmax, alpha)

    print('drdv0 ~', ((r2-r1)/(2*delta) - drdv0.dot(dvec)) / drdv0.dot(dvec))
    print('dvdv0 ~', ((v2-v1)/(2*delta) - dvdv0.dot(dvec)) / dvdv0.dot(dvec))
    print('dtdv0 ~', ((t2-t1)/(2*delta) - dtdv0.dot(dvec)) / dtdv0.dot(dvec))

    r1,v1,t1,_ = fandg(r0, v0, p - dvec * delta, tau, kmax, alpha)
    r2,v2,t2,_ = fandg(r0, v0, p + dvec * delta, tau, kmax, alpha)

    print('drdp  ~', ((r2-r1)/(2*delta) - drdp.dot(dvec)) / drdp.dot(dvec))
    print('dvdp  ~', ((v2-v1)/(2*delta) - dvdp.dot(dvec)) / dvdp.dot(dvec))
    print('dtdp  ~', ((t2-t1)/(2*delta) - dtdp.dot(dvec)) / dtdp.dot(dvec))

if __name__ == "__main__":
    main()
