:
:  ichanWT2005.mod 
:
:   Alan Goldin Lab, University of California, Irvine
:   Jay Lickfett - Last Modified: 6 July 2005
:
:  This file is the Nav1.1 wild-type channel model described in:
:
:		Barela et al. An Epilepsy Mutation in the Sodium Channel SCN1A That Decreases
:	    Channel Excitability.  J. Neurosci. 26(10): p. 2714-2723 
:
:
:   The model is derived from the one described in:
:
:    	Spampanato et al. (2004a) Increased Neuronal Firing in Computer Simulations 
:		of Sodium Channel Mutations that Cause Generalized Epilepsy with Febrile Seizures Plus.
:		Journal of Neurophysiology 91:2040-2050
:
:	and
:
:	 	Spampanato et al. (2004b) A Novel Epilepsy Mutation 
:   	in the Sodium Channel SCN1A Identifies a Cytoplasmic Domain for 
:		Beta Subunit Interaction. J. Neurosci. 24(44):10022-10034
:
 

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (uF) = (microfarad)
    (molar) = (1/liter)
    (nA) = (nanoamp)
    (mM) = (millimolar)
    (um) = (micron)
    (S) = (siemens)
    FARADAY = 96520 (coul)
    R = 8.3134  (joule/degC)

}

NEURON { 
    SUFFIX inak2005 
    USEION nat READ enat WRITE inat VALENCE 1
    USEION kf READ ekf WRITE ikf  VALENCE 1
:    NONSPECIFIC_CURRENT il 
    RANGE gnat, gkf
    RANGE gnatbar, gkfbar, gnablock
:    RANGE gl, el f
    RANGE minf, mtau, hshift, sshift, mvhalf, mk, hvhalf, hk, svhalf, sk, mtaubase, htauk, htauvhalf, htauk, stauvhalf, stauk, hinf, htau, sinf, stau, nfinf, nftau, inat, m, h, s, htaubase, staubase
}

 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

 
PARAMETER {
    
    celsius = 6.3 (degC)
    dt (ms) 
    enat  (mV)
    gnatbar (mho/cm2)   
    ekf  (mV)
    gkfbar (mho/cm2)
    type = 0       : 0 is WT, 1 is T875M, 2 is W1204R, 3 is R1648H, 4 is R859C
    gnablock = 1.0
    mvhalf = 27.4 (mV)
    mk = 5.4043
    hvhalf = 41.9 (mV)
    hk = 6.7
    svhalf = 46.0 (mV)
    sk = 6.6
    hshift = 0 (mV)
    sshift = 0 (mV)
    mtaubase = 0.15 (ms)
    htaubase = 23.12 (ms)
    htauvhalf = 77.58 (mV)
    htauk = 43.92
    staubase = 140400 (ms)
    stauvhalf = 71.3 (mV)
    stauk = 30.9
}


ASSIGNED {
      
    v (mV) 
    gnat (mho/cm2) 
    gkf (mho/cm2)
    inat (mA/cm2)
    ikf (mA/cm2)
    minf hinf sinf nfinf
    mtau (ms) htau (ms) stau (ms) nftau (ms)
    mexp hexp sexp nfexp
} 


STATE {
    m h s nf
}
 

BREAKPOINT {
    SOLVE states
    gnat = gnatbar*gnablock*m*m*m*h*s  
    inat = gnat*(v - enat)
    gkf = gkfbar*nf*nf*nf*nf
    ikf = gkf*(v-ekf)
}

 
UNITSOFF

 
INITIAL {

    trates(v)
    
    m = minf
    h = hinf
    s = sinf

    nf = nfinf
}


PROCEDURE states() {        : Computes state variables m, h, s and n 
                            : at the current v and dt.        
    trates(v)           

    m = m + mexp*(minf-m)
    h = h + hexp*(hinf-h)
    s = s + sexp*(sinf-s)
    nf = nf + nfexp*(nfinf-nf)
}
 

LOCAL q10


PROCEDURE rates(v (mV)) {   :Computes rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, sum, vhs, vss
    q10 = 3^((celsius - 6.3)/10)
    
    vhs = v - hshift
    vss = v - sshift

    :"m" sodium activation system
    minf = 1/(1+exp(-(v+mvhalf)/mk))		
    mtau = mtaubase

    :"h" sodium fast inactivation system
    hinf = 1/(1+exp((vhs+hvhalf)/hk))				
    htau = htaubase*exp(-0.5*((v+htauvhalf)/htauk)^2) 	
       
    :"s" sodium slow inactivation system
    sinf = 1/(1+exp((vss+svhalf)/sk))						
    stau = staubase*exp(-0.5*((v+stauvhalf)/stauk)^2)	

    :"nf" fKDR activation system				
    alpha = -0.07*vtrap((v+65-47),-6)
    beta = 0.264/exp((v+65-22)/40)
    sum = alpha+beta        
    nftau = 1/sum      
    nfinf = alpha/sum
}
 

PROCEDURE trates(v (mV)) {  :Build table with rate and other constants at current v.
                            :Call once from HOC to initialize inf at resting v.
    LOCAL tinc
 
    TABLE minf, mexp, hinf, hexp, sinf, sexp, nfinf, nfexp, mtau, htau, stau, nftau
        DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
    rates(v)    : not consistently executed from here if usetable_hh == 1
                : so don't expect the tau values to be tracking along with
                : the inf values in hoc

    tinc = -dt * q10
    mexp = 1 - exp(tinc/mtau)
    hexp = 1 - exp(tinc/htau)
    sexp = 1 - exp(tinc/stau)
    nfexp = 1 - exp(tinc/nftau)
}

 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.

    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{  
        vtrap = x/(exp(x/y) - 1)
    }
}
 

UNITSON

