
"""
Rank-2 anisotropic rotational diffusion propagator.

This module implements the closed-form rank-2 Wigner-basis propagator for
rotational diffusion with principal diffusion constants Dx, Dy, Dz.

Basis ordering:
    m = [2, 1, 0, -1, -2]

The propagated matrix E(t) is defined by

    <D^{2*}_{a b}(0) D^2_{c d}(t)> = (1/5) * delta_{a c} * E_{b d}(t)

for an isotropic equilibrium ensemble.

The code automatically switches to exact isotropic and axially symmetric
limits when the corresponding diffusion constants are equal to within a
tolerance. This avoids divisions by zero and minimizes cancellation errors.
"""


import math
import numpy as np
import sys

M_ORDER = (2, 1, 0, -1, -2)


def _sinhc(x: float) -> float:
    """Return sinh(x)/x with a stable small-x expansion."""
    ax = abs(x)
    if ax < 1e-8:
        x2 = x * x
        # 1 + x^2/6 + x^4/120 + x^6/5040
        return 1.0 + x2 / 6.0 + (x2 * x2) / 120.0 + (x2 * x2 * x2) / 5040.0
    return math.sinh(x) / x


def _is_close(a: float, b: float, tol: float) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= tol * scale



def IsIsotropic(
    Dx: float,
    Dy: float,
    Dz: float,
    tol: float = 1e-12
    )-> bool:
    # Exact isotropic limit: Dx = Dy = Dz
    if _is_close(Dx, Dy, tol) and _is_close(Dx, Dz, tol):
        return True
    return False

def IsAxial(
    Dx: float,
    Dy: float,
    Dz: float,
    tol: float = 1e-12
    )-> bool:
    # Exact axial symmetry: Dx = Dy (about z)
    if _is_close(Dx, Dy, tol):
        return True
    return False



def rank2_generator_matrix(Dx: float, Dy: float, Dz: float) -> np.ndarray:
    """
    Return the 5x5 rank-2 generator H in the basis [2, 1, 0, -1, -2].

    The propagator is E(t) = exp(-H t).
    """
    Dx = float(Dx)
    Dy = float(Dy)
    Dz = float(Dz)

    S = Dx + Dy
    delta = Dx - Dy
    A = S + 4.0 * Dz
    B = 2.5 * S + Dz
    C = 3.0 * S
    U = (math.sqrt(6.0) / 2.0) * delta
    V = 1.5 * delta

    H = np.array(
        [
            [A, 0.0, U, 0.0, 0.0],
            [0.0, B, 0.0, V, 0.0],
            [U, 0.0, C, 0.0, U],
            [0.0, V, 0.0, B, 0.0],
            [0.0, 0.0, U, 0.0, A],
        ],
        dtype=float,
    )
    return H



def rank2_propagator(
    t: float,
    Dx: float,
    Dy: float,
    Dz: float,
    tol: float = 1e-12,
) -> np.ndarray:
    """
    Closed-form 5x5 propagator E(t) = exp(-H t) for rank 2.

    Parameters
    ----------
    t : float
        Time.
    Dx, Dy, Dz : float
        Principal-axis rotational diffusion coefficients.
    tol : float
        Relative tolerance for deciding isotropic or axially symmetric limits.

    Returns
    -------
    E : (5, 5) ndarray
        Propagator in the basis [2, 1, 0, -1, -2].
    """
    t = float(t)
    Dx = float(Dx)
    Dy = float(Dy)
    Dz = float(Dz)

    if t == 0.0:
        return np.eye(5, dtype=float)

    # Exact isotropic limit: Dx = Dy = Dz
    #if _is_close(Dx, Dy, tol) and _is_close(Dx, Dz, tol):
    if(IsIsotropic(Dx,Dy,Dz,tol)):
        return np.exp(-6.0 * Dx * t) * np.eye(5, dtype=float)

    # Exact axial symmetry: Dx = Dy (about z)
    #if _is_close(Dx, Dy, tol):
    if(IsAxial(Dx,Dy,Dz,tol)):
        Dperp = 0.5 * (Dx + Dy)
        Dpar = Dz
        E = np.zeros((5, 5), dtype=float)

        # Basis order: m = [2, 1, 0, -1, -2]
        e2 = math.exp(-(2.0 * Dperp + 4.0 * Dpar) * t)
        e1 = math.exp(-(5.0 * Dperp + Dpar) * t)
        e0 = math.exp(-6.0 * Dperp * t)

        E[0, 0] = e2
        E[1, 1] = e1
        E[2, 2] = e0
        E[3, 3] = e1
        E[4, 4] = e2
        return E

    # Fully anisotropic closed form.
    S = Dx + Dy
    delta = Dx - Dy
    A = S + 4.0 * Dz
    B = 2.5 * S + Dz
    mu = 2.0 * (S + Dz)

    alpha = 2.0 * Dz - S
    beta = math.sqrt(3.0) * delta
    kappa = math.sqrt(alpha * alpha + beta * beta)
    omega = 1.5 * delta

    # Stable versions of the potentially singular expressions:
    #   sinh(kappa t)/kappa = t * sinhc(kappa t)
    shc = _sinhc(kappa * t)
    k_sinh = t * shc

    expA = math.exp(-A * t)
    expB = math.exp(-B * t)
    expMu = math.exp(-mu * t)

    coshK = math.cosh(kappa * t)
    sinhW = math.sinh(omega * t)
    coshW = math.cosh(omega * t)

    alpha_term = alpha * k_sinh
    beta_term = beta * k_sinh / math.sqrt(2.0)

    E = np.zeros((5, 5), dtype=float)

    # Even block: m = 2, 0, -2 mixed by anisotropy
    x = expMu * (coshK - alpha_term)
    y = expMu * (coshK + alpha_term)

    E[0, 0] = 0.5 * (expA + x)
    E[4, 4] = E[0, 0]

    E[0, 4] = 0.5 * (x - expA)
    E[4, 0] = E[0, 4]

    E[0, 2] = -beta_term * expMu
    E[2, 0] = E[0, 2]
    E[4, 2] = E[0, 2]
    E[2, 4] = E[0, 2]

    E[2, 2] = y

    # Odd block: m = 1, -1
    E[1, 1] = expB * coshW
    E[3, 3] = E[1, 1]
    E[1, 3] = -expB * sinhW
    E[3, 1] = E[1, 3]


    
    return E


def Jomega(T,om):
    return T/(1+T**2*om**2)


def rank2_propagatorDict(
    Dx: float,
    Dy: float,
    Dz: float,
    tol: float = 1e-12,
) -> (dict,dict):
    """
    Closed-form 5x5 propagator E(t) = exp(-H t) for rank 2.
    Returns diction EE with refactors and TT with timescales.
    Both are index -2 to +2 but including only non-zero elements.
    The propgator E(t) is recovered by taking EE*exp(-TT*t), looping
    over all elements. The timescales need to be kept separate as
    for NMR they will find themselves in effective timescales
    in spectral density functions.

    The setup is clever enough to detect if we are isotropic, axially
    symmetric or fully anisotropic.

    For isotropic or axial, the dictionary EE contains only diagonal elements
    with a single timescale for each (as Wigners are the eigenfunctions).
    For fully anisotropic this is no longer the case, and there is more
    mixing.

    For axial, requires Dz be the unique axis, Dx/Dy are degenerate.

    EE and TT are recombined into E(t) using function AssembleE
    
    Parameters
    ----------
    Dx, Dy, Dz : float
        Principal-axis rotational diffusion coefficients.
    tol : float
        Relative tolerance for deciding isotropic or axially symmetric limits.

    Returns
    -------
    EE,TT : indexed [i][j], both -2 to +2 but containing non-zero elements. 
    """
    #t = float(t)
    Dx = float(Dx)
    Dy = float(Dy)
    Dz = float(Dz)

    #blank the dictionaries
    EE={};TT={}
    for i in range(5):
        EE[i]={}
        TT[i]={}

    # Exact isotropic limit: Dx = Dy = Dz
    #if _is_close(Dx, Dy, tol) and _is_close(Dx, Dz, tol):
    if(IsIsotropic(Dx,Dy,Dz,tol)):
        for i in range(5):
            EE[i][i]=(1,);
            TT[i][i]=(6.0 * Dx ,)
        return EE,TT   #completed isotropic case.
        #E=np.eye(5,dtype=float)
        #T=E*6.0 * Dx 
        #return np.exp(-6.0 * Dx * t) * np.eye(5, dtype=float)
        #return E,T

    # Exact axial symmetry: Dx = Dy (about z)
    #if _is_close(Dx, Dy, tol):
    if(IsAxial(Dx,Dy,Dz,tol)):
        Dperp = 0.5 * (Dx + Dy)
        Dpar = Dz
        E = np.zeros((5, 5), dtype=float)
        T = np.zeros((5, 5), dtype=float)

        # Basis order: m = [2, 1, 0, -1, -2]
        e2 = 1 #math.exp(-(2.0 * Dperp + 4.0 * Dpar) * t)
        e1 = 1 #math.exp(-(5.0 * Dperp + Dpar) * t)
        e0 = 1 #math.exp(-6.0 * Dperp * t)

        for i in range(5):
            EE[i][i]=(1,);

        TT[0][0] = (2.0 * Dperp + 4.0 * Dpar,)
        TT[1][1] = (5.0 * Dperp + Dpar,)
        TT[2][2] = (6.0 * Dperp,)
        TT[3][3] = (5.0 * Dperp + Dpar,)
        TT[4][4] = (2.0 * Dperp + 4.0 * Dpar,)

        return EE,TT

    #def _sinhc(x: float) -> float:
    #    """Return sinh(x)/x with a stable small-x expansion."""
    #    ax = abs(x)
    #    if ax < 1e-8:
    #        x2 = x * x
    #        # 1 + x^2/6 + x^4/120 + x^6/5040
    #        return 1.0 + x2 / 6.0 + (x2 * x2) / 120.0 + (x2 * x2 * x2) / 5040.0
    #    return math.sinh(x) / x
    #shc = _sinhc(kappa * t)
    #shc = 0.5*( math.exp( kappa * t) -  math.exp( -kappa * t) )
    #k_sinh = t * shc
    #we are not doing the small value expansion here. Seems to be no need? 
    
    
    # Fully anisotropic closed form.
    S = Dx + Dy
    delta = Dx - Dy
    A = S + 4.0 * Dz
    B = 2.5 * S + Dz
    mu = 2.0 * (S + Dz)

    alpha = 2.0 * Dz - S
    beta = math.sqrt(3.0) * delta
    kappa = math.sqrt(alpha * alpha + beta * beta)
    omega = 1.5 * delta

    
    #keep the dictionary for the prefactors EE and the rate constants TT seperate
    EE[0][0]=[];TT[0][0]=[]
    EE[0][0].append(0.5);                   TT[0][0].append(A)
    EE[0][0].append(0.25*(1-alpha/kappa));  TT[0][0].append(mu-kappa)
    EE[0][0].append(0.25*(1+alpha/kappa));  TT[0][0].append(mu+kappa)

    EE[4][4]=EE[0][0]
    TT[4][4]=TT[0][0]

    EE[0][4]=[];TT[0][4]=[]
    EE[0][4].append(-0.5);                 TT[0][4].append(A)
    EE[0][4].append(0.25*(1-alpha/kappa)); TT[0][4].append(mu-kappa)
    EE[0][4].append(0.25*(1+alpha/kappa)); TT[0][4].append(mu+kappa)

    EE[4][0]=EE[0][4]
    TT[4][0]=TT[0][4]

    EE[0][2]=[];TT[0][2]=[]
    EE[0][2].append( (-beta/kappa  * 0.25 *math.sqrt(2.0)) );    TT[0][2].append(mu-kappa)
    EE[0][2].append( (+beta/kappa  * 0.25 *math.sqrt(2.0)) );    TT[0][2].append(mu+kappa)

    EE[2][0]=EE[0][2];  TT[2][0]=TT[0][2]
    EE[4][2]=EE[0][2];  TT[4][2]=TT[0][2]
    EE[2][4]=EE[0][2];  TT[2][4]=TT[0][2]

    EE[2][2]=[];TT[2][2]=[]
    EE[2][2].append( 0.5*(1+alpha/kappa) );    TT[2][2].append(mu-kappa)
    EE[2][2].append( 0.5*(1-alpha/kappa) );    TT[2][2].append(mu+kappa)

    EE[1][1]=[];TT[1][1]=[]
    EE[1][1].append( 0.5 );    TT[1][1].append(B-omega)
    EE[1][1].append( 0.5 );    TT[1][1].append(B+omega)

    EE[3][3]=EE[1][1]; TT[3][3]=TT[1][1]

    EE[1][3]=[];TT[1][3]=[]
    EE[1][3].append( -0.5 );    TT[1][3].append(B-omega)
    EE[1][3].append(  0.5 );    TT[1][3].append(B+omega)
    
    EE[3][1]=EE[1][3]; TT[3][1]=TT[1][3]

    return EE,TT  #return the two dictionaries.
    

    """
    #Write E as a sum of exponentials to note prefactors and timescales
    E = np.zeros((5, 5), dtype=float)
    E[0, 0] = 0.5 * math.exp(-A * t)    +     0.25* math.exp( (-mu +kappa) * t) + 0.25*math.exp( (-mu -kappa) * t)  -    0.25*alpha/kappa  * math.exp( (-mu +kappa) * t)     + 0.25*alpha/kappa * math.exp( (-mu -kappa) * t)      
    E[4, 4] = E[0, 0]
    E[0, 4] =  0.25*math.exp((-mu+kappa) * t)   + 0.25* math.exp((-mu-kappa) * t) -    alpha/kappa  * 0.25* math.exp((-mu+kappa) * t) +  alpha/kappa  * 0.25* math.exp((-mu-kappa) * t)    -    0.5 *math.exp(-A * t) 
    E[4, 0] = E[0, 4]
    E[0, 2] =  -beta/kappa  * 0.5 / math.sqrt(2.0) * math.exp( (-mu + kappa) * t) + beta/kappa  * 0.5 / math.sqrt(2.0) * math.exp((-mu-kappa) * t)   
    E[2, 0] = E[0, 2]
    E[4, 2] = E[0, 2]
    E[2, 4] = E[0, 2]
    E[2, 2] =  0.5*math.exp((-mu+kappa) * t)   +    0.5*math.exp((-mu-kappa) * t)   +  alpha/kappa  * 0.5* math.exp(( -mu+ kappa) * t) - alpha/kappa  * 0.5* math.exp( ( -mu-kappa) * t) 
    # Odd block: m = 1, -1
    E[1, 1] =   0.5*math.exp( (-B+omega) * t) + 0.5*math.exp((-B-omega) * t)
    E[3, 3] = E[1, 1]
    E[1, 3] = -  0.5* math.exp((-B+omega) * t) +0.5*math.exp((-B-omega) * t)
    E[3, 1] = E[1, 3]
    return E
    """
    



#take trace of matrix. Not sure why, chatGPT liked this.
def rank2_scalar_correlation(
    t: float,
    Dx: float,
    Dy: float,
    Dz: float,
    tol: float = 1e-12,
) -> float:
    """
    Scalar trace-average correlation C(t) = (1/5) Tr[E(t)].
    """
    E = rank2_propagator(t, Dx, Dy, Dz, tol=tol)
    return float(np.trace(E) / 5.0)


def rank2_autocorrelation_tensor_element(
    t: float,
    a: int,
    b: int,
    c: int,
    d: int,
    Dx: float,
    Dy: float,
    Dz: float,
    tol: float = 1e-12,
) -> float:
    """
    Return <D^{2*}_{a b}(0) D^2_{c d}(t)> for the rank-2 model.

    Indices a,c are the space-fixed indices and b,d are the body-fixed indices.
    No cacheing is done: this is very brute force. For testing only.
    
    """
    if a != c:
        return 0.0

    try:
        ib = M_ORDER.index(b)
        id_ = M_ORDER.index(d)
    except ValueError as exc:
        raise ValueError(f"b and d must be in {M_ORDER}") from exc

    E = rank2_propagator(t, Dx, Dy, Dz, tol=tol)
    return float(E[ib, id_] / 5.0)

def rank2_autocorrelation_tensor_elementStr(
    a: int,
    b: int,
    c: int,
    d: int,
    FULLY: bool
) -> str:
    """
    Return a string for b and d elements that will later be numerically
    evalulated to give a symbolic element for b/d coefficients.

    Later on the symoblic elements will be turned into <D^{2*}_{a b}(0) D^2_{c d}(t)> for the rank-2 model.

    the A and C elements provide a delta function, which will be incorporated
    explicitly.

    It returns false if the cominbation of b and d will be zero.
    
    Indices a,c are the space-fixed indices and b,d are the body-fixed indices.

    FULLY indicates if we are fully anisotropic (True) or if we are isotropic or axial (False)
    
    """
    if a != c:
        return False

    try:
        ib = M_ORDER.index(b)
        id_ = M_ORDER.index(d)
    except ValueError as exc:
        raise ValueError(f"b and d must be in {M_ORDER}") from exc

    #E = rank2_propagator(t, Dx, Dy, Dz, tol=tol)
    #return float(E[ib, id_] / 5.0)

    #add some more check to verify that we don't have a zero value
    if(ib==id_): #this is always good.
        return "h("+str(ib)+","+str(id_)+")"

    if(FULLY==False):
        return False
    good=True
    #if fully anisotropic, then there are options.
    if(ib==0 or ib==2 or ib==4):
        if id_ not in (0,2,4):
            good=False
    if(ib==1 or ib==3):
        if id_ not in (1,3):
            good=False
    if(good==False):
        return False
    #print(ib,id_,good)
    
    return "h("+str(ib)+","+str(id_)+")"


def calcRank2fromStr(
        h: str,
        E: np.ndarray
) -> float:
    """
    Return <D^{2*}_{a b}(0) D^2_{c d}(t)> for the rank-2 model.
    Take matrix E, and return just the right elements. Unneccessary, really.
    But useful for testing.
    
    Indices a,c are the space-fixed indices and b,d are the body-fixed indices.
    """
    if(h==False):
        return 0.0
    
    ib=int(h.split('(')[1].split(',')[0])
    id_=int(h.split(',')[1].split(')')[0])

    """
    try:
        ib = M_ORDER.index(b)
        id_ = M_ORDER.index(d)
    except ValueError as exc:
        raise ValueError(f"b and d must be in {M_ORDER}") from exc
    """
    #E = rank2_propagator(t, Dx, Dy, Dz, tol=tol)

    #return rank2_propagatorTime(t, Dx, Dy, Dz, tol=tol)
    
    return float(E[ib, id_] / 5.0)


def AssembleE(
        t: float,
        EE: dict,
        TT: dict
              )->float:
    #reassemble the correlation function dictionary
    E = np.zeros((5, 5), dtype=float)
    for key,vals in EE.items():
        for koi,vols in vals.items():
            for V,T in zip(vols,TT[key][koi]):
                E[key,koi]+=V*math.exp(-T*t)
    return E

def CalcJ(
        om: float,
        EE: dict,
        TT: dict
)->float:
    #compute spectral density function
    J=0.0
    for key,vals in EE.items():
        for koi,vols in vals.items():
            for V,T in zip(vols,TT[key][koi]):
                J+=V* Jomega(T,om)
    #print("J:",J)
    return J



def validate_rank2_implementation() -> None:
    """
    Run numerical checks.

    Checks:
      1. isotropic limit
      2. axially symmetric limit
      3. comparison against expm(-H t) if SciPy is available
    """
    import numpy.testing as npt

    # 1) Isotropic check
    t = 0.37
    D = 1.23
    E_iso = rank2_propagator(t, D, D, D)
    E_iso_expected = math.exp(-6.0 * D * t) * np.eye(5)
    npt.assert_allclose(E_iso, E_iso_expected, rtol=1e-14, atol=1e-14)

    # 2) Axially symmetric check
    Dperp = 0.8
    Dpar = 1.9
    E_ax = rank2_propagator(t, Dperp, Dperp, Dpar)
    E_ax_expected = np.diag(
        [
            math.exp(-(2.0 * Dperp + 4.0 * Dpar) * t),
            math.exp(-(5.0 * Dperp + Dpar) * t),
            math.exp(-6.0 * Dperp * t),
            math.exp(-(5.0 * Dperp + Dpar) * t),
            math.exp(-(2.0 * Dperp + 4.0 * Dpar) * t),
        ]
    )
    npt.assert_allclose(E_ax, E_ax_expected, rtol=1e-14, atol=1e-14)

    # 3) General case check against expm, if SciPy is available
    try:
        from scipy.linalg import expm
    except Exception:
        return

    Dx, Dy, Dz = 1.4, 0.6, 2.2
    H = rank2_generator_matrix(Dx, Dy, Dz)
    E_ref = expm(-H * t)
    E_closed = rank2_propagator(t, Dx, Dy, Dz)
    npt.assert_allclose(E_closed, E_ref, rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    validate_rank2_implementation()
    print("All rank-2 rotational diffusion checks passed.")

    #from rank2_anisotropic_rotational_diffusion import rank2_propagator, rank2_scalar_correlation

    Dx=1.2
    Dy=0.8
    #Dy=1.2
    Dz=2.0
    #Dz=1.2
    
    E = rank2_propagator(t=0.3, Dx=Dx, Dy=Dy, Dz=Dz)         #compute E matrix numerically
    print(E)


    
    #C = rank2_scalar_correlation(t=0.3, Dx=Dx, Dy=Dy, Dz=Dz) 
    #print(C)


    #need to separate prefactors and time scales for relaxation calcs.
    EE,TT=rank2_propagatorDict(Dx=Dx, Dy=Dy, Dz=Dz)  #get numerical dictionary of prefactors (EE) and timescales (TT)

    ET=AssembleE(t=0.3,EE=EE,TT=TT)  #reassemble E matrix from dictionary at specified time (for testing)
    

    residual=np.sum(np.fabs(E-ET)) #work out difference between numerical, and reassembled E matrices
    print('residual:',residual)    #note: in relcalc we will never use this, but we will use EE and TT.
    if(residual>1E-9):             #      we check here that EE and TT are correct.
        print('shit')
        print(E)
        print(ET)
        sys.exit(100)


    FULLY=False    #set to true for fully anisotropic. set to false for either isotropic or axially symmetric.
                   #this preconsiders which terms in <D^*_ab D_cd> are zero. For False, only b=d terms. For true, more mixing is allowed.
                   # note: a=c is ALWAYS true. The difference between fully anistropic and not is that more b/d mixing is allowed.
                   # In recalc, we will assign a string to the expressions that links to the b/d values as a place holder
                   # When we calculate numerical rates, we will numerically evaluate these factors via these protocols.
                   # The recipie will be: first calculate EE and TT, dictionaries for prefactors and constants that depend on DxDyDz
                   # Then for each h(b,d) in a relaxation rate, just the relevant row/column.
                   # The similarity of Dx/Dy and Dy/Dz are checked so that if we say we are fully anisoptric, we relax nicely
                   # to the appropriate limit.

                   
    ISO=IsIsotropic(Dx,Dy,Dz)
    AX=IsAxial(Dx,Dy,Dz)
    print("Isotropic?",ISO)
    print("Axial?",AX)
    if(FULLY==False and ISO==False and AX==False):
        print("Error: we have fully anisotropic values, but are not using the fully anisotropic model")
        print("Changing FULLY to true")
        FULLY=True  #change and crack on.

    if(FULLY==True and (ISO==True or AX==True)):
        print("Note we have isotropic values, but we are still using fully anisotropic model")
        print("This is fine,  it will collapse correctly, but you are wasting FLOPS!")
        FULLY=False #manually over-ride

         
                   
    a=1;c=1;   #set integers for a and c.
    
    #first print the allowed strings.
    for b in (-2,-1,0,1,2):
        for d in (-2,-1,0,1,2):
            eleStr=rank2_autocorrelation_tensor_elementStr(a=1,b=b,c=1,d=d,FULLY=FULLY)
            if(eleStr!=False):
                print(eleStr)

    #now, for given and c, calculate the E matrix. verify symbolic and numerical treatments are the same.
    #in RelCalc, we will retain the symbolc factors from rank2_autocorrelation_tensor_elementStr
    #and turn them into numerical values with refactors/timecale as part of the numerical engine.
    for b in (-2,-1,0,1,2):
        for d in (-2,-1,0,1,2):
            print()
            eleStr=rank2_autocorrelation_tensor_elementStr(a=1,b=b,c=1,d=d,FULLY=FULLY)
            print('string element:',eleStr,'(',b,',',d,')')

            ele=rank2_autocorrelation_tensor_element(t=0.3,a=1,b=b,c=1,d=d,Dx=Dx,Dy=Dy,Dz=Dz)
            eleStrVal=calcRank2fromStr(eleStr,ET)

            print(ele,eleStrVal)
            if(math.fabs(eleStrVal-ele)>1E-8):  #check residual for errors.
                print(b,d,'shit')
                sys.exit(100)
            
            

            J=CalcJ(om=4.0,EE=EE,TT=TT)
            print("J:",J)
            
            #eleStrVal





