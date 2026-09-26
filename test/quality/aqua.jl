using Aqua
using GeometricIntegratorsBase
using Test

Aqua.test_all(
    GeometricIntegratorsBase;
    undefined_exports = (broken = true,),  # issue #39: initialguess! is exported, not defined
    stale_deps = false,                    # issue #40: run as @test_broken below
    piracies = (broken = true,)            # issue #41: integrate and integrate! on foreign types
)

# `test_stale_deps` takes no `broken` keyword, so its check is run here directly.
# issue #40: SafeTestsets is in [deps] of Project.toml and nothing under src/ loads it
@testset "Stale dependencies" begin
    @test_broken isempty(Aqua.find_stale_deps(Base.PkgId(GeometricIntegratorsBase)))  # issue #40
end
