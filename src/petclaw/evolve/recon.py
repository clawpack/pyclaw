#Reconstruction functions for SharpClaw
import numpy as np

def weno5(q):

    epweno=1.e-36

    dqiph = np.diff(q,1,0)

    LL=3
    UL=q.shape[0]-2
    qr=np.empty(q.shape)
    ql=np.empty(q.shape)

    for m1 in [1,2]:
        #m1=1: construct q^-_{i+1/2} (ql)
        #m1=2: construct q^+_{i+1/2} (qr)
        im=(-1)**(m1+1)
        ione=im
        inone=-im
        intwo=-2*im

        #Create references to DQ slices
        dq_intwo=dqiph[LL+intwo-1:UL+intwo-1,:]
        dq_ione =dqiph[LL+ione-1 :UL+ione-1 ,:]
        dq_inone=dqiph[LL+inone-1:UL+inone-1,:]
        dq      =dqiph[LL-1:UL-1            ,:]

        t1 = im*(dq_intwo-dq_inone)
        t2 = im*(dq_inone-dq)
        t3 = im*(dq      -dq_ione)

        tt1=13.*t1**2+3.*(   dq_intwo - 3.*dq_inone)**2
        tt2=13.*t2**2+3.*(   dq_inone +    dq      )**2
        tt3=13.*t3**2+3.*(3.*dq       -    dq_ione )**2

        tt1=(epweno+tt1)**2
        tt2=(epweno+tt2)**2
        tt3=(epweno+tt3)**2
        s1 = tt2*tt3
        s2 = 6.*tt1*tt3
        s3 = 3.*tt1*tt2
        t0 = 1./(s1+s2+s3)
        s1 *= t0
        s3 *= t0

        z=(s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3. \
                + (-q[LL-2:UL-2,:]+7.*(q[LL-1:UL-1,:]+q[LL:UL,:])-q[LL+1:UL+1,:])/12.
        if m1==1: qr[LL-1:UL-1,:] = z
        else: ql[LL:UL,:] = z
                
    return ql,qr

def weno5_wave(q,wave,s):

    epweno=1.e-36

    qr=q.copy()
    ql=q.copy()
    LL=2
    UL=q.shape[0]-3
    mwaves=wave.shape[2]
    meqn=wave.shape[1]
    for m1 in [1,2]:
        #m1=1: construct q^-_{i+1/2} (ql)
        #m1=2: construct q^+_{i+1/2} (qr)
        im=(-1)**(m1+1)
        ione=im
        inone=-im
        intwo=-2*im

        for mw in xrange(mwaves):
            wnorm2 = wave[LL:UL,0,mw]**2
            theta1 = wave[LL+intwo:UL+intwo,0,mw]*wave[LL:UL,0,mw]
            theta2 = wave[LL+inone:UL+inone,0,mw]*wave[LL:UL,0,mw]
            theta3 = wave[LL+ione :UL+ione ,0,mw]*wave[LL:UL,0,mw]
            for m in xrange(1,meqn):
                wnorm2 += wave[LL:UL,m,mw]**2
                theta1 += wave[LL+intwo:UL+intwo,m,mw]*wave[LL:UL,m,mw]
                theta2 += wave[LL+inone:UL+inone,m,mw]*wave[LL:UL,m,mw]
                theta3 += wave[LL+ione :UL+ione ,m,mw]*wave[LL:UL,m,mw]

            t1=im*(theta1-theta2)
            t2=im*(theta2-wnorm2)
            t3=im*(wnorm2-theta3)

            tt1=13.*t1**2+3.*(   theta1 - 3.*theta2)**2
            tt2=13.*t2**2+3.*(   theta2 +    wnorm2)**2
            tt3=13.*t3**2+3.*(3.*wnorm2 -    theta3)**2

            tt1=(epweno+tt1)**2
            tt2=(epweno+tt2)**2
            tt3=(epweno+tt3)**2
            s1 = tt2*tt3
            s2 = 6.*tt1*tt3
            s3 = 3.*tt1*tt2
            t0 = 1./(s1+s2+s3)
            s1 *= t0
            s3 *= t0

            z=(s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3. \
                    + im*(theta2+6.*wnorm2-theta3)/12.
            u=np.where(wnorm2>1.e-14,z,0.)
            wnorm2=np.where(wnorm2>1.e-14,1./wnorm2,1.)

            for m in xrange(meqn):
                if m1==1: qr[LL:UL,m] += u*wave[LL:UL,m,mw]*wnorm2
                else: ql[LL+1:UL+1,m] += u*wave[LL:UL,m,mw]*wnorm2
                
    return ql,qr
