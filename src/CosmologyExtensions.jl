"""
    pk_from_sigma8(z, h, Ωb, Ωm, ΩΛ, ns, σ8, Tcmb)
    returns a tidy function of comoving wavenumber in units of inverse Mpc
                        to return normalised linear P(k) in units of Mpc^3
    based on astro-ph/9710252 and astro-ph/9709112
    Tcmb in units of K
"""
const G(y) = y*(-6*sqrt(1+y)+(2+3*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)))
const amin = 1e-3
const W = kR->3*(sin(kR)/kR^3-cos(kR)/kR^2)
const Wsq = kR->9*(sin(kR)/kR^3-cos(kR)/kR^2)^2
const dWdkR = kR->9*cos(kR)/kR^3-(9/kR^2-3)*sin(kR)/kR^2

function σ2_quadgk(pk, R::Float64)
    QuadGK.quadgk(k->Wsq(k*R)*pk(k)*k^2/2/pi^2,0,20/R)[1]
end

function dσ2dR_quadgk(pk, R::Float64)
    QuadGK.quadgk(k->W(k*R)*dWdkR(k*R)*pk(k)*k^3/pi^2,0,20/R)[1]
end

function tinker10bM(lnν::Float64,δc::Float64=1.686,Δm::Float64=200.)
    y = log10(Δm)
    A = 1 + 0.24*y*exp(-256/y^4)
    a = 0.44*y - 0.88
    B = 0.183 # fit
    b = 1.5 # fit
    C = 0.019 + 0.107*y + 0.19*exp(-256/y^4)
    c = 2.4 # fit
    return 1 - A*exp(lnν*a)/(exp(lnν*a)+δc^a) + B*exp(lnν*b) + C*exp(lnν*c)
end

function tEH_setup(z::Float64,h::Float64,Ωb::Float64,Ωm::Float64,ΩΛ::Float64,ns::Float64,σ8::Float64,Tcmb::Float64)
    fb = Ωb/Ωm
    Ω0h2 = Ωm*h^2
    Ωbh2 = Ωb*h^2
    Θcmb = Tcmb/2.7
    b1d = 0.313*Ω0h2^(-0.419)*(1+0.607*Ω0h2^0.674)
    b2d = 0.238*Ω0h2^0.223
    z_d = 1291*(Ω0h2)^(0.251)/(1+0.659*Ω0h2^0.828)*(1+b1d*Ωbh2^b2d)
    # note on zeq:
    # the preprint clearly asserts that the quantity below is zeq
    # but equally clearly, the actual code available from background.uchicago.edu asserts it is 1+zeq
    # therefore https://github.com/steven-murray/hmf goes with the code
    # we're going with the code too
    zeqp1 = 2.5e4*Ω0h2*Θcmb^(-4)
    keq = 7.46e-2*Ω0h2*Θcmb^(-2) # in inverse Mpc, still no h
    # OK, another note, this time on this baryon-photon ratio R
    # once again the preprint clearly asserts that the correct dependence is (z/1e3)^{-1}
    # once again http://background.uchicago.edu/~whu/transfer/power.f clearly asserts it is 1e3/(1+z)
    # once again hmf goes with the code
    # so once again we're going with the code
    R_d = 31.5*Ωbh2*Θcmb^(-4)/((z_d+1)/1e3)
    Req = 31.5*Ωbh2*Θcmb^(-4)/(zeqp1/1e3)
    # note that since both these redshift values are quite high, the difference in practice is a few percent at most
    # but by god it matters when you're trying to validate code

    s = 2/(3*keq)*sqrt(6/Req)*log((sqrt(1+R_d)+sqrt(R_d+Req))/(1+sqrt(Req)))

    a1 = (46.9*Ω0h2)^0.670*(1+(32.1*Ω0h2)^(-0.532))
    a2 = (12.0*Ω0h2)^0.424*(1+(45.0*Ω0h2)^(-0.582))
    αc = a1^(-fb)*a2^(-fb^3)
    b1 = 0.944/(1+(458*Ω0h2)^(-0.708))
    b2 = (0.395*Ω0h2)^(-0.0266)
    βc = 1/(1+b1*((1-fb)^b2-1))

    kSilk = 1.6*Ωbh2^0.52*Ω0h2^0.73*(1+(10.4*Ω0h2)^(-0.95)) # again in inverse Mpc without h

    q = Θcmb^2/(Ω0h2)
    αb = 2.07*keq*s*(1+R_d)^(-3/4)*G(zeqp1/(1+z_d))
    βb = 0.5+fb+(3-2*fb)*sqrt((17.2*Ω0h2)^2+1)
    βnode = 8.41*Ω0h2^0.435
    a = 1/(1+z)
    fb, Ω0h2, Ωbh2, Θcmb, b1d, b2d, z_d, R_d, zeqp1, keq, Req, s, a1, a2, αc, b1, b2, βc, kSilk, q, αb, βb, βnode, a
end

# a quick note on g^2(z):
# the original formula from astro-ph/9710252 looks like this:
# gsq = a->ΩΛ+Ωm/a^3+(1-Ωm-ΩΛ)/a^2
# but we recognise that this basically assumes everything that isn't DE or matter goes to curvature
# which ... isn't the best assumption??? you're way better off assuming the rest is radiation or massless neutrinos
# so we will use this:
const gsq = (a,ΩΛ,Ωm)->ΩΛ+Ωm/a^3+(1-Ωm-ΩΛ)/a^4 

function pk_from_sigma8(z::Float64,h::Float64,Ωb::Float64,Ωm::Float64,ΩΛ::Float64,ns::Float64,σ8::Float64,Tcmb::Float64)
    fb, Ω0h2, Ωbh2, Θcmb, b1d, b2d, z_d, R_d, zeqp1, keq, Req, s, a1, a2, αc, b1, b2, βc, kSilk, q, αb, βb, βnode, a = tEH_setup(z,h,Ωb,Ωm,ΩΛ,ns,σ8,Tcmb)#::NTuple{24,Float64}
    T0tilde_1 = k->log(ℯ+1.8*βc*(q*k))/(log(ℯ+1.8*βc*(q*k))+(14.2+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    T0tilde_αc = k->log(ℯ+1.8*βc*(q*k))/(log(ℯ+1.8*βc*(q*k))+(14.2/αc+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    Tc = k->T0tilde_1(k)/(1+(k*s/5.4)^4)+(1-1/(1+(k*s/5.4)^4))*T0tilde_αc(k)

    T0tilde_11 = k->log(ℯ+1.8*(q*k))/(log(ℯ+1.8*(q*k))+(14.2+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    Tb = k->(T0tilde_11(k)/(1+(k*s/5.2)^2)+αb/(1+(βb/k/s)^3)*exp(-(k/kSilk)^1.4))*sin(k*s/(1+(βnode/k/s)^3)^(1/3))/(k*s/(1+(βnode/k/s)^3)^(1/3))

    D1r = (sqrt(gsq(a,ΩΛ,Ωm))*QuadGK.quadgk(a->1/(a^3*gsq(a,ΩΛ,Ωm)^(3/2)),amin,a)[1]/QuadGK.quadgk(a->1/(a^3*gsq(a,ΩΛ,Ωm)^(3/2)),amin,1)[1])
    # calculate sigma8 to find normalisation at z = 0; note that the EH transfer function is z-independent
    σ8sq_wo_norm = σ2_quadgk(k->k^ns*(fb*Tb(k) + (1-fb)*Tc(k))^2,8/h)
    totnorm = (D1r*σ8)^2/σ8sq_wo_norm
    k -> totnorm*k^ns*(fb*Tb(k) + (1-fb)*Tc(k))^2
end

function pk_and_hmf_from_sigma8(z::Float64,h::Float64,Ωb::Float64,Ωm::Float64,ΩΛ::Float64,ns::Float64,σ8::Float64,Tcmb::Float64,R::Vector{Float64})
    fb, Ω0h2, Ωbh2, Θcmb, b1d, b2d, z_d, R_d, zeqp1, keq, Req, s, a1, a2, αc, b1, b2, βc, kSilk, q, αb, βb, βnode, a = tEH_setup(z,h,Ωb,Ωm,ΩΛ,ns,σ8,Tcmb)#::NTuple{24,Float64}
    T0tilde_1 = k->log(ℯ+1.8*βc*(q*k))/(log(ℯ+1.8*βc*(q*k))+(14.2+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    T0tilde_αc = k->log(ℯ+1.8*βc*(q*k))/(log(ℯ+1.8*βc*(q*k))+(14.2/αc+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    Tc = k->T0tilde_1(k)/(1+(k*s/5.4)^4)+(1-1/(1+(k*s/5.4)^4))*T0tilde_αc(k)

    T0tilde_11 = k->log(ℯ+1.8*(q*k))/(log(ℯ+1.8*(q*k))+(14.2+386/(1+69.9*(q*k)^1.08))*(q*q*k*k))
    Tb = k->(T0tilde_11(k)/(1+(k*s/5.2)^2)+αb/(1+(βb/k/s)^3)*exp(-(k/kSilk)^1.4))*sin(k*s/(1+(βnode/k/s)^3)^(1/3))/(k*s/(1+(βnode/k/s)^3)^(1/3))

    D1r = (sqrt(gsq(a,ΩΛ,Ωm))*QuadGK.quadgk(a->1/(a^3*gsq(a,ΩΛ,Ωm)^(3/2)),amin,a)[1]/QuadGK.quadgk(a->1/(a^3*gsq(a,ΩΛ,Ωm)^(3/2)),amin,1)[1])
    # calculate sigma8 to find normalisation at z = 0; note that the EH transfer function is z-independent
    σ8sq_wo_norm = σ2_quadgk(k->k^ns*(fb*Tb(k) + (1-fb)*Tc(k))^2,8/h)
    totnorm = (D1r*σ8)^2/σ8sq_wo_norm
    # and now the hmf part
    nR = size(R)[1]
    hmf_σ2 = Vector{Float64}(undef,nR)
    hmf_dlnσ2dlnR = Vector{Float64}(undef,nR)
    hmf_dndlnM = Vector{Float64}(undef,nR)
    
    for i = 1:nR
        hmf_σ2[i] = quadde(k->Wsq(k*R[i])*totnorm*k^(2+ns)*(fb*Tb(k) + (1-fb)*Tc(k))^2/2/pi^2,0,20/R[i])[1]
        hmf_dlnσ2dlnR[i] = R[i] / hmf_σ2[i] * quadde(k->W(k*R[i])*dWdkR(k*R[i])*totnorm*k^(3+ns)*(fb*Tb(k) + (1-fb)*Tc(k))^2/pi^2,0,20/R[i])[1]
        hmf_dndlnM[i] = -hmf_dlnσ2dlnR[i] * HaloMF.tinker08MF(2*log(1.686) - log(hmf_σ2[i]), z, 200) / (4*pi*R[i]^3) # (Mh_range[i]*3/ρm)
    end
    k -> totnorm*k^ns*(fb*Tb(k) + (1-fb)*Tc(k))^2, hmf_σ2, hmf_dlnσ2dlnR, log(1.686) .- log.(hmf_σ2)./2., hmf_dndlnM
end

function pk_and_hmf_from_sigma8!(z::Float64,h::Float64,Ωb::Float64,Ωm::Float64,ΩΛ::Float64,ns::Float64,σ8::Float64,Tcmb::Float64,R::Vector{Float64},hmf_σ2::Vector{Float64},lnν::Vector{Float64},hmf_dlnσ2dlnR::Vector{Float64},hmf_dndlnM::Vector{Float64},hmf_dndM::Vector{Float64})
    fb, Ω0h2, Ωbh2, Θcmb, b1d, b2d, z_d, R_d, zeqp1, keq, Req, s, a1, a2, αc, b1, b2, βc, kSilk, q, αb, βb, βnode, a = tEH_setup(z,h,Ωb,Ωm,ΩΛ,ns,σ8,Tcmb)#::NTuple{24,Float64}

    T0tilde = (q,α,β)->log(ℯ+1.8*β*(q))/(log(ℯ+1.8*β*(q))+(14.2/α+386.0/(1+69.9*(q)^1.08))*(q*q))
    Tb = k->(T0tilde(q*k,1,1)/(1+(k*s/5.2)^2)+αb/(1+(βb/k/s)^3)*exp(-(k/kSilk)^1.4))*sin(k*s/(1+(βnode/k/s)^3)^(1/3))/(k*s/(1+(βnode/k/s)^3)^(1/3))
    Tc = k->T0tilde(q*k,1,βc)/(1+(k*s/5.4)^4)+(1-1/(1+(k*s/5.4)^4))*T0tilde(q*k,αc,βc)

    D1r = (sqrt(gsq(a,ΩΛ,Ωm))*quadde(a->1/(a^3*gsq(a,ΩΛ,Ωm)^(3/2)),amin,a)[1]/quadde(a->1/(a^3*gsq(a,ΩΛ,Ωm)^(3/2)),amin,1)[1])
    # calculate sigma8 to find normalisation at z = 0; note that the EH transfer function is z-independent
    σ8sq_wo_norm = σ2_quadgk(k->k^ns*(fb*Tb(k) + (1-fb)*Tc(k))^2,8.0/h)
    totnorm = (D1r*σ8)^2/σ8sq_wo_norm
    # and now the hmf part
    nR = size(R)[1]
    for i = 1:nR
        hmf_σ2[i] = quadde(k->Wsq(k*R[i])*totnorm*k^(2+ns)*(fb*Tb(k) + (1-fb)*Tc(k))^2/2.0/pi^2,0.0,20.0/R[i])[1]
        lnν[i] = log(1.686) - log(hmf_σ2[i])/2
        hmf_dlnσ2dlnR[i] = R[i] / hmf_σ2[i] * quadde(k->W(k*R[i])*dWdkR(k*R[i])*totnorm*k^(3+ns)*(fb*Tb(k) + (1-fb)*Tc(k))^2/pi^2,0.0,20.0/R[i])[1]
        hmf_dndlnM[i] = -hmf_dlnσ2dlnR[i] * HaloMF.tinker08MF(2*lnν[i], z, 200) / (4*pi*R[i]^3) # (Mh_range[i]*3/ρm)
        hmf_dndM[i] = hmf_dndlnM[i] * 3 / (Ω0h2*2.775e11*4*pi*R[i]^3)
    end
    k -> totnorm*k^ns*(fb*Tb(k) + (1-fb)*Tc(k))^2
end
