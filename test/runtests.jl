using LIMahl
using Test

@testset "LIMahl.jl" begin
    using Cosmology, Unitful, UnitfulAstro
    pkp_cosmo = Cosmology.cosmology(h=0.7,OmegaM=0.286,Tcmb=2.725)
    testmodel = LIMahl.LineModel(115.0*Unitful.GHz,30.0*Unitful.GHz,pkp_cosmo,0.96,0.82,0.047,2.725,1e9*UnitfulAstro.Msun,1e15*UnitfulAstro.Msun,200,M->2e-6*UnitfulAstro.Lsun*(M/UnitfulAstro.Msun));
    # normalisation can be kind of shaky so suppose we just want within 3% of the lim result
    @test isapprox(testmodel.Pk(0.01035), 8.93e3, rtol=0.03)
    @test isapprox(testmodel.Tmean, 15.283 * Unitful.μK, rtol=0.03)
    @test isapprox(testmodel.Pshot, 17539. * Unitful.μK^2*UnitfulAstro.Mpc^3, rtol=0.03)
end
