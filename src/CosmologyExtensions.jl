"""
    pk_from_sigma8(z, h, Ωb, Ωm, ΩΛ, ns, σ8, Tcmb)
    returns a tidy function of comoving wavenumber in units of inverse Mpc
                        to return normalised linear P(k) in units of Mpc^3
    based on astro-ph/9710252 and astro-ph/9709112
    Tcmb in units of K
"""
G(y) = y*(-6*sqrt(1+y)+(2+3*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)))
amin = 1e-3

function pk_from_sigma8(z::Real,h::Real,Ωb::Real,Ωm::Real,ΩΛ::Real,ns::Real,σ8::Real,Tcmb::Real)
    gsq = a->ΩΛ+Ωm/a^3+(1-Ωm-ΩΛ)/a^2
    fb = Ωb/Ωm
    Ω0h2 = Ωm*h^2
    Ωbh2 = Ωb*h^2
    Θcmb = Tcmb/2.7
    b1d = 0.313*Ω0h2^(-0.419)*(1+0.607*Ω0h2^0.674)
    b2d = 0.238*Ω0h2^0.223
    z_d = 1291*(Ω0h2)^(0.251)/(1+0.659*Ω0h2^0.828)*(1+b1d*Ωbh2^b2d)
    R_d = 31.5*Ωbh2*Θcmb^(-4)/(z_d/1e3)
    zeq = 2.5e4*Ω0h2*Θcmb^(-4)
    keq = 7.48e-2*Ω0h2*Θcmb^(-2) # in inverse Mpc, still no h
    Req = 31.5*Ωbh2*Θcmb^(-4)/(zeq/1e3)

    s = 2/(3*keq)*sqrt(6/Req)*log((sqrt(1+R_d)+sqrt(R_d+Req))/(1+sqrt(Req)))

    a1 = (46.9*Ω0h2)^0.670*(1+(32.1*Ω0h2)^(-0.532))
    a2 = (12.0*Ω0h2)^0.424*(1+(45.0*Ω0h2)^(-0.582))
    αc = a1^(-fb)*a2^(-fb^3)
    b1 = 0.944/(1+(458*Ω0h2)^(-0.708))
    b2 = (0.395*Ω0h2)^(-0.0266)
    βc = 1/(1+b1*((1-fb)^b2-1))

    kSilk = 1.6*Ωbh2^0.52*Ω0h2^0.73*(1+(10.4*Ω0h2)^(-0.95)) # again in inverse Mpc without h

    q = (Tcmb/2.7)^2/(Ω0h2)
    T0tilde_1 = k->log(ℯ+1.8*βc*(q*k))/(log(ℯ+1.8*βc*(q*k))+(14.2+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    T0tilde_αc = k->log(ℯ+1.8*βc*(q*k))/(log(ℯ+1.8*βc*(q*k))+(14.2/αc+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    Tc = k->T0tilde_1(k)/(1+(k*s/5.4)^4)+(1-1/(1+(k*s/5.4)^4))*T0tilde_αc(k)

    T0tilde_11 = k->log(ℯ+1.8*(q*k))/(log(ℯ+1.8*(q*k))+(14.2+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    αb = 2.07*keq*s*(1+R_d)^(-3/4)*G((1+zeq)/(1+z_d))
    βb = 0.5+fb+(3-2*fb)*sqrt((17.2*Ω0h2)^2+1)
    βnode = 8.41*Ω0h2^0.435
    Tb = k->(T0tilde_11(k)/(1+(k*s/5.2)^2)+αb/(1+(βb/k/s)^3)*exp(-(k/kSilk)^1.4))*sin(k*s/(1+(βnode/k/s)^3)^(1/3))/(k*s/(1+(βnode/k/s)^3)^(1/3))

    a = 1/(1+z)
    D1r = sqrt(gsq(1/(1+z)))*QuadGK.quadgk(a->1/(a^3*gsq(a)^(3/2)),amin,1/(1+z))[1]/QuadGK.quadgk(a->1/(a^3*gsq(a)^(3/2)),amin,1)[1]
    # calculate sigma8 to find normalisation at z = 0; note that the EH transfer function is z-independent
    σ8sq_wo_norm = MatterPower.sigma2(k->k^ns*(fb*Tb(k) + (1-fb)*Tc(k))^2,8/h)
    totnorm = (D1r*σ8)^2/σ8sq_wo_norm
    k -> totnorm*k^ns*(fb*Tb(k) + (1-fb)*Tc(k))^2
end
