# Create a Bouncing Ball FMU

Tutorial by Johannes Stoljar, Tobias Thummerer, Simon Exner | Last edit: October 29 2024

## License


```julia
# Copyright (c) 2021 Tobias Thummerer, Lars Mikelsons, Josef Kircher, Johannes Stoljar
# Licensed under the MIT license.
# See LICENSE (https://github.com/thummeto/FMIExport.jl/blob/main/LICENSE) file in the project root for details.
```

## Motivation

This Julia Package FMIExport.jl is motivated by the export of simulation models in Julia. Here the FMI specification is implemented. FMI (Functional Mock-up Interface) is a free standard ([fmi-standard.org](https://fmi-standard.org)) that defines a container and an interface to exchange dynamic models using a combination of XML files, binaries and C code zipped into a single file. The user is able to create own FMUs (Functional Mock-up Units).

## Target group

The example is primarily intended for users who work in the field of simulations. The example wants to show how simple it is to export FMUs in Julia.

## Introduction to the example

This example shows how to export a FMU from julia-code. It uses the BouncingBall FMU, that can be found on the main branch of FMIExport in [examples/FMI2/BouncingBall](https://github.com/ThummeTo/FMIExport.jl/tree/main/examples/FMI2/BouncingBall). This notebook will show you how to export it.

## Installation prerequisites

|     | Description                       | Command                   | Alternative                                    |
|:----|:----------------------------------|:--------------------------|:-----------------------------------------------|
| 1.  | Enter Package Manager via         | ]                         |                                                |
| 2.  | Install FMIExport via             | add FMIExport             | add "https://github.com/ThummeTo/FMIExport.jl" |
| 3.  | Install FMIBuild via              | add FMIBuild              | add "https://github.com/ThummeTo/FMIBuild.jl"  |

## REPL-commands or build-script

The way to do this usually will be the REPL, but if you plan on exporting FMUs in an automated way, you may want to use a jl script containing the following commands.
To run this example, the previously installed packages must be included.


```julia
using FMIExport
using FMIBuild: saveFMU
```

next we have to define where to put the generated files


```julia
tmpDir = mktempdir(; prefix="fmibuildjl_test_", cleanup=false) 
@info "Saving example files at: $(tmpDir)"
fmu_save_path = joinpath(tmpDir, "BouncingBall.fmu")  
```

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39mSaving example files at: C:\Users\RUNNER~1\AppData\Local\Temp\fmibuildjl_test_TfcMf5
    




    "C:\\Users\\RUNNER~1\\AppData\\Local\\Temp\\fmibuildjl_test_TfcMf5\\BouncingBall.fmu"



Remember, that we use the FMU-source stored at [examples/FMI2/BouncingBall](https://github.com/ThummeTo/FMIExport.jl/tree/main/examples/FMI2/BouncingBall). If you execute this notebook locally, make shure to ajust the fmu_source_path to where your FMU-Package resides. **It is important, that an absolute path is provided!** For this notebook to work in the automated bulid pipeline, this absolute path is obtained by the following instructions. If you run this example locally, you can provide the path manually, just make shure you use the correct directory seperator or just use just use julias `joinpath` function.


```julia
working_dir = pwd() # current working directory
println(string("pwd() returns: ", working_dir))

package_dir = split(working_dir, joinpath("examples", "jupyter-src"))[1] # remove everything after and including "examples\jupyter-src"
println(string("package_dir is ", package_dir))

fmu_source_package = joinpath(package_dir, "examples", "FMI2", "BouncingBall") # add correct relative path
println(string("fmu_source_package is ", fmu_source_package))

fmu_source_path = joinpath(fmu_source_package, "src", "BouncingBall.jl") # add correct relative path
println(string("fmu_source_path is ", fmu_source_path))
```

    pwd() returns: D:\a\FMIExport.jl\FMIExport.jl\examples\jupyter-src
    package_dir is D:\a\FMIExport.jl\FMIExport.jl\
    fmu_source_package is D:\a\FMIExport.jl\FMIExport.jl\examples\FMI2\BouncingBall
    fmu_source_path is D:\a\FMIExport.jl\FMIExport.jl\examples\FMI2\BouncingBall\src\BouncingBall.jl
    

TODO The following codecell contains *workardound* code that will be obsolete with the next release


```julia
using FMIExport.FMIBase.FMICore: fmi2True, fmi2False 

EPS = 1e-6

FMU_FCT_INIT = function()
    m = 1.0         # ball mass
    r = 0.0         # ball radius
    d = 0.7         # ball collision damping
    v_min = 1e-1    # ball minimum velocity
    g = 9.81        # gravity constant 
    sticking = fmi2False

    s = 1.0         # ball position
    v = 0.0         # ball velocity
    a = 0.0         # ball acceleration

    t = 0.0        
    x_c = [s, v]      
    ẋ_c = [v, a]
    x_d = [sticking]
    u = []
    p = [m, r, d, v_min, g]

    return (t, x_c, ẋ_c, x_d, u, p)
end

FMU_FCT_EVALUATE = function(t, x_c, ẋ_c, x_d, u, p, eventMode)
    m, r, d, v_min, g = p
    s, v = x_c
    sticking = x_d[1]
    _, a = ẋ_c

    if sticking == fmi2True
        a = 0.0
    elseif sticking == fmi2False
        if eventMode
            if s < r && v < 0.0
                s = r + EPS # so that indicator is not triggered again
                v = -v*d 
                
                # stop bouncing to prevent high frequency bouncing (and maybe tunneling the floor)
                if abs(v) < v_min
                    sticking = fmi2True
                    v = 0.0
                end
            end
        else
            # no specials in continuos time mode
        end

        a = (m * -g) / m     # the system's physical equation (a little longer than necessary)
    else
        @error "Unknown value for `sticking` == $(sticking)."
        return (x_c, ẋ_c, x_d, p)
    end

    x_c = [s, v]
    ẋ_c = [v, a]
    x_d = [sticking]
    p = [m, r, d, v_min, g]

    return (x_c, ẋ_c, x_d, p) # evaluation can't change discrete state!
end

FMU_FCT_OUTPUT = function(t, x_c, ẋ_c, x_d, u, p)
    m, r, d, v_min, g = p
    s, v = x_c
    _, a = ẋ_c
    sticking = x_d[1]

    y = [s]

    return y
end

FMU_FCT_EVENT = function(t, x_c, ẋ_c, x_d, u, p)
    m, r, d, v_min, g = p
    s, v = x_c
    _, a = ẋ_c
    sticking = x_d[1]
   
    if sticking == fmi2True
        z1 = 1.0            # event 1: ball stay-on-ground
    else
        z1 = (s-r)          # event 1: ball hits ground 
    end

    z = [z1]

    return z
end
FMIBUILD_CONSTRUCTOR = function(resPath="")
    fmu = fmi2CreateSimple(initializationFct=FMU_FCT_INIT,
                        evaluationFct=FMU_FCT_EVALUATE,
                        outputFct=FMU_FCT_OUTPUT,
                        eventFct=FMU_FCT_EVENT)

    fmu.modelDescription.modelName = "BouncingBall"

    # modes 
    fmi2ModelDescriptionAddModelExchange(fmu.modelDescription, "BouncingBall")

    # states [2]
    fmi2AddStateAndDerivative(fmu, "ball.s"; stateDescr="Absolute position of ball center of mass", derivativeDescr="Absolute velocity of ball center of mass")
    fmi2AddStateAndDerivative(fmu, "ball.v"; stateDescr="Absolute velocity of ball center of mass", derivativeDescr="Absolute acceleration of ball center of mass")

    # discrete state [1]
    fmi2AddIntegerDiscreteState(fmu, "sticking"; description="Indicator (boolean) if the mass is sticking on the ground, as soon as abs(v) < v_min")

    # outputs [1]
    fmi2AddRealOutput(fmu, "ball.s_out"; description="Absolute position of ball center of mass")

    # parameters [5]
    fmi2AddRealParameter(fmu, "m";     description="Mass of ball")
    fmi2AddRealParameter(fmu, "r";     description="Radius of ball")
    fmi2AddRealParameter(fmu, "d";     description="Collision damping constant (velocity fraction after hitting the ground)")
    fmi2AddRealParameter(fmu, "v_min"; description="Minimal ball velocity to enter on-ground-state")
    fmi2AddRealParameter(fmu, "g";     description="Gravity constant")

    fmi2AddEventIndicator(fmu)

    return fmu
end
fmu = FMIBUILD_CONSTRUCTOR()
```




    Model name:	BouncingBall
    Type:		0



TODO? It is questionable if this is the job of the library or the user... Currently it is not implemented and therefor the job of the user

We need to make shure the fmu_source_package is instantiated


```julia
using Pkg
notebook_env = Base.active_project(); # save current enviroment to return to it after we are done
Pkg.activate(fmu_source_package); # activate the FMUs enviroment

# make shure to use the same FMI source as in the enviroment of this example ("notebook_env"). 
# As this example is automattically built using the local FMIExport package and not the one from the Juila registry, we need to add it using "develop". 
Pkg.develop(PackageSpec(path=package_dir)); # If you added FMIExport using "add FMIExport", you have to remove this line and use instantiate instead.
# Pkg.instantiate(); # instantiate the FMUs enviroment only if develop was not previously called

Pkg.activate(notebook_env); # return to the original notebooks enviroment
```

    [32m[1m  Activating[22m[39m project at `D:\a\FMIExport.jl\FMIExport.jl\examples\FMI2\BouncingBall`

    
    

    [32m[1m   Resolving[22m[39m package versions...
    

    [32m[1m    Updating[22m[39m `D:\a\FMIExport.jl\FMIExport.jl\examples\FMI2\BouncingBall\Project.toml`
      [90m[226f0e26] [39m[92m+ FMIBuild v0.3.2[39m
      [90m[31b88311] [39m[92m+ FMIExport v0.4.0 `D:\a\FMIExport.jl\FMIExport.jl\`[39m
    [32m[1m    Updating[22m[39m `D:\a\FMIExport.jl\FMIExport.jl\examples\FMI2\BouncingBall\Manifest.toml`
    

      [90m[47edcb42] [39m[92m+ ADTypes v1.9.0[39m
      [90m[7d9f7c33] [39m[92m+ Accessors v0.1.38[39m
      [90m[79e6a3ab] [39m[92m+ Adapt v4.1.1[39m
      [90m[4fba245c] [39m[92m+ ArrayInterface v7.16.0[39m
      [90m[4c555306] [39m[92m+ ArrayLayouts v1.10.4[39m
      [90m[62783981] [39m[92m+ BitTwiddlingConvenienceFunctions v0.1.6[39m
      [90m[2a0fbf3d] [39m[92m+ CPUSummary v0.2.6[39m
    [33m⌅[39m [90m[d360d2e6] [39m[92m+ ChainRulesCore v1.24.0[39m
      [90m[fb6a15b2] [39m[92m+ CloseOpenIntervals v0.1.13[39m
      [90m[38540f10] [39m[92m+ CommonSolve v0.2.4[39m
      [90m[bbf7d656] [39m[92m+ CommonSubexpressions v0.3.1[39m
      [90m[f70d9fcc] [39m[92m+ CommonWorldInvalidations v1.0.0[39m
      [90m[34da2185] [39m[92m+ Compat v4.16.0[39m
      [90m[a33af91c] [39m[92m+ CompositionsBase v0.1.2[39m
      [90m[2569d6c7] [39m[92m+ ConcreteStructs v0.2.3[39m
      [90m[187b0558] [39m[92m+ ConstructionBase v1.5.8[39m
      [90m[adafc99b] [39m[92m+ CpuId v0.3.1[39m
      [90m[9a962f9c] [39m[92m+ DataAPI v1.16.0[39m
      [90m[864edb3b] [39m[92m+ DataStructures v0.18.20[39m
      [90m[e2d170a0] [39m[92m+ DataValueInterfaces v1.0.0[39m
      [90m[2b5f629d] [39m[92m+ DiffEqBase v6.158.3[39m
    [33m⌅[39m [90m[459566f4] [39m[92m+ DiffEqCallbacks v3.9.1[39m
      [90m[163ba53b] [39m[92m+ DiffResults v1.1.0[39m
      [90m[b552c78f] [39m[92m+ DiffRules v1.15.1[39m
      [90m[a0c0ee7d] [39m[92m+ DifferentiationInterface v0.6.18[39m
      [90m[ffbed154] [39m[92m+ DocStringExtensions v0.9.3[39m
      [90m[4e289a0a] [39m[92m+ EnumX v1.0.4[39m
      [90m[f151be2c] [39m[92m+ EnzymeCore v0.8.5[39m
      [90m[e2ba6199] [39m[92m+ ExprTools v0.1.10[39m
    [33m⌅[39m [90m[6b7a57c9] [39m[92m+ Expronicon v0.8.5[39m
      [90m[8f5d6c58] [39m[92m+ EzXML v1.2.0[39m
      [90m[900ee838] [39m[92m+ FMIBase v1.0.10[39m
      [90m[226f0e26] [39m[92m+ FMIBuild v0.3.2[39m
      [90m[8af89139] [39m[92m+ FMICore v1.1.1[39m
      [90m[31b88311] [39m[92m+ FMIExport v0.4.0 `D:\a\FMIExport.jl\FMIExport.jl\`[39m
      [90m[7034ab61] [39m[92m+ FastBroadcast v0.3.5[39m
      [90m[9aa1b823] [39m[92m+ FastClosures v0.3.2[39m
      [90m[29a986be] [39m[92m+ FastLapackInterface v2.0.4[39m
      [90m[1a297f60] [39m[92m+ FillArrays v1.13.0[39m
      [90m[6a86dc24] [39m[92m+ FiniteDiff v2.26.0[39m
      [90m[f6369f11] [39m[92m+ ForwardDiff v0.10.37[39m
      [90m[069b7b12] [39m[92m+ FunctionWrappers v1.1.3[39m
      [90m[77dc65aa] [39m[92m+ FunctionWrappersWrappers v0.1.3[39m
      [90m[d9f16b24] [39m[92m+ Functors v0.4.12[39m
    [32m⌃[39m [90m[46192b85] [39m[92m+ GPUArraysCore v0.1.6[39m
      [90m[c27321d9] [39m[92m+ Glob v1.3.1[39m
      [90m[3e5b6fbb] [39m[92m+ HostCPUFeatures v0.1.17[39m
      [90m[615f187c] [39m[92m+ IfElse v0.1.1[39m
      [90m[3587e190] [39m[92m+ InverseFunctions v0.1.17[39m
      [90m[92d709cd] [39m[92m+ IrrationalConstants v0.2.2[39m
      [90m[82899510] [39m[92m+ IteratorInterfaceExtensions v1.0.0[39m
      [90m[692b3bcd] [39m[92m+ JLLWrappers v1.6.1[39m
      [90m[ef3ab10e] [39m[92m+ KLU v0.6.0[39m
      [90m[ba0b0d4f] [39m[92m+ Krylov v0.9.8[39m
      [90m[10f19ff3] [39m[92m+ LayoutPointers v0.1.17[39m
      [90m[5078a376] [39m[92m+ LazyArrays v2.2.1[39m
      [90m[87fe0de2] [39m[92m+ LineSearch v0.1.4[39m
      [90m[d3d80556] [39m[92m+ LineSearches v7.3.0[39m
      [90m[7ed4a6bd] [39m[92m+ LinearSolve v2.36.2[39m
      [90m[2ab3a3ac] [39m[92m+ LogExpFunctions v0.3.28[39m
      [90m[bdcacae8] [39m[92m+ LoopVectorization v0.12.171[39m
      [90m[d8e11817] [39m[92m+ MLStyle v0.4.17[39m
      [90m[1914dd2f] [39m[92m+ MacroTools v0.5.13[39m
      [90m[d125e4d3] [39m[92m+ ManualMemory v0.1.8[39m
      [90m[bb5d69b7] [39m[92m+ MaybeInplace v0.1.4[39m
      [90m[46d2c3a1] [39m[92m+ MuladdMacro v0.2.4[39m
      [90m[d41bc354] [39m[92m+ NLSolversBase v7.8.3[39m
      [90m[77ba4419] [39m[92m+ NaNMath v1.0.2[39m
    [33m⌅[39m [90m[8913a72c] [39m[92m+ NonlinearSolve v3.15.1[39m
      [90m[6fe1bfb0] [39m[92m+ OffsetArrays v1.14.1[39m
      [90m[bac558e1] [39m[92m+ OrderedCollections v1.6.3[39m
      [90m[9b87118b] [39m[92m+ PackageCompiler v2.1.22[39m
      [90m[65ce6f38] [39m[92m+ PackageExtensionCompat v1.0.2[39m
      [90m[d96e819e] [39m[92m+ Parameters v0.12.3[39m
      [90m[f517fe37] [39m[92m+ Polyester v0.7.16[39m
      [90m[1d0040c9] [39m[92m+ PolyesterWeave v0.2.2[39m
      [90m[d236fae5] [39m[92m+ PreallocationTools v0.4.24[39m
      [90m[aea7be01] [39m[92m+ PrecompileTools v1.2.1[39m
      [90m[21216c6a] [39m[92m+ Preferences v1.4.3[39m
      [90m[92933f4c] [39m[92m+ ProgressMeter v1.10.2[39m
      [90m[3cdcf5f2] [39m[92m+ RecipesBase v1.3.4[39m
      [90m[731186ca] [39m[92m+ RecursiveArrayTools v3.27.2[39m
      [90m[f2c3362d] [39m[92m+ RecursiveFactorization v0.2.23[39m
      [90m[189a3867] [39m[92m+ Reexport v1.2.2[39m
      [90m[05181044] [39m[92m+ RelocatableFolders v1.0.1[39m
      [90m[ae029012] [39m[92m+ Requires v1.3.0[39m
      [90m[7e49a35a] [39m[92m+ RuntimeGeneratedFunctions v0.5.13[39m
      [90m[94e857df] [39m[92m+ SIMDTypes v0.1.0[39m
      [90m[476501e8] [39m[92m+ SLEEFPirates v0.6.43[39m
      [90m[0bca4576] [39m[92m+ SciMLBase v2.58.0[39m
      [90m[19f34311] [39m[92m+ SciMLJacobianOperators v0.1.1[39m
      [90m[c0aeaf25] [39m[92m+ SciMLOperators v0.3.12[39m
      [90m[53ae85a6] [39m[92m+ SciMLStructures v1.5.0[39m
      [90m[6c6a2e73] [39m[92m+ Scratch v1.2.1[39m
      [90m[efcf1570] [39m[92m+ Setfield v1.1.1[39m
    [33m⌅[39m [90m[727e6d20] [39m[92m+ SimpleNonlinearSolve v1.12.3[39m
      [90m[9f842d2f] [39m[92m+ SparseConnectivityTracer v0.6.8[39m
      [90m[0a514795] [39m[92m+ SparseMatrixColorings v0.4.8[39m
      [90m[e56a9233] [39m[92m+ Sparspak v0.3.9[39m
      [90m[276daf66] [39m[92m+ SpecialFunctions v2.4.0[39m
      [90m[aedffcd0] [39m[92m+ Static v1.1.1[39m
      [90m[0d7ed370] [39m[92m+ StaticArrayInterface v1.8.0[39m
      [90m[1e83bf80] [39m[92m+ StaticArraysCore v1.4.3[39m
      [90m[7792a7ef] [39m[92m+ StrideArraysCore v0.5.7[39m
      [90m[2efcf032] [39m[92m+ SymbolicIndexingInterface v0.3.34[39m
      [90m[3783bdb8] [39m[92m+ TableTraits v1.0.1[39m
      [90m[bd369af6] [39m[92m+ Tables v1.12.0[39m
      [90m[8290d209] [39m[92m+ ThreadingUtilities v0.5.2[39m
      [90m[a759f4b9] [39m[92m+ TimerOutputs v0.5.25[39m
      [90m[d5829a12] [39m[92m+ TriangularSolve v0.2.1[39m
      [90m[781d530d] [39m[92m+ TruncatedStacktraces v1.4.0[39m
      [90m[3a884ed6] [39m[92m+ UnPack v1.0.2[39m
      [90m[3d5dd08c] [39m[92m+ VectorizationBase v0.21.71[39m
      [90m[a5390f91] [39m[92m+ ZipFile v0.10.1[39m
      [90m[1d5cc7b8] [39m[92m+ IntelOpenMP_jll v2024.2.1+0[39m
      [90m[94ce4f54] [39m[92m+ Libiconv_jll v1.17.0+1[39m
      [90m[856f044c] [39m[92m+ MKL_jll v2024.2.0+0[39m
      [90m[efe28fd5] [39m[92m+ OpenSpecFun_jll v0.5.5+0[39m
      [90m[02c8fc9c] [39m[92m+ XML2_jll v2.13.4+0[39m
      [90m[1317d2d5] [39m[92m+ oneTBB_jll v2021.12.0+0[39m
      [90m[0dad84c5] [39m[92m+ ArgTools v1.1.1[39m
      [90m[56f22d72] [39m[92m+ Artifacts[39m
      [90m[2a0f44e3] [39m[92m+ Base64[39m
      [90m[ade2ca70] [39m[92m+ Dates[39m
      [90m[8ba89e20] [39m[92m+ Distributed[39m
      [90m[f43a241f] [39m[92m+ Downloads v1.6.0[39m
      [90m[7b1f6079] [39m[92m+ FileWatching[39m
      [90m[9fa8497b] [39m[92m+ Future[39m
      [90m[b77e0a4c] [39m[92m+ InteractiveUtils[39m
      [90m[4af54fe1] [39m[92m+ LazyArtifacts[39m
      [90m[b27032c2] [39m[92m+ LibCURL v0.6.4[39m
      [90m[76f85450] [39m[92m+ LibGit2[39m
      [90m[8f399da3] [39m[92m+ Libdl[39m
      [90m[37e2e46d] [39m[92m+ LinearAlgebra[39m
      [90m[56ddb016] [39m[92m+ Logging[39m
      [90m[d6f4376e] [39m[92m+ Markdown[39m
      [90m[ca575930] [39m[92m+ NetworkOptions v1.2.0[39m
      [90m[44cfe95a] [39m[92m+ Pkg v1.10.0[39m
      [90m[de0858da] [39m[92m+ Printf[39m
      [90m[3fa0cd96] [39m[92m+ REPL[39m
      [90m[9a3f8284] [39m[92m+ Random[39m
      [90m[ea8e919c] [39m[92m+ SHA v0.7.0[39m
      [90m[9e88b42a] [39m[92m+ Serialization[39m
      [90m[6462fe0b] [39m[92m+ Sockets[39m
      [90m[2f01184e] [39m[92m+ SparseArrays v1.10.0[39m
      [90m[10745b16] [39m[92m+ Statistics v1.10.0[39m
      [90m[fa267f1f] [39m[92m+ TOML v1.0.3[39m
      [90m[a4e569a6] [39m[92m+ Tar v1.10.0[39m
      [90m[8dfed614] [39m[92m+ Test[39m
      [90m[cf7118a7] [39m[92m+ UUIDs[39m
      [90m[4ec0a83e] [39m[92m+ Unicode[39m
      [90m[e66e0078] [39m[92m+ CompilerSupportLibraries_jll v1.1.1+0[39m
      [90m[deac9b47] [39m[92m+ LibCURL_jll v8.4.0+0[39m
      [90m[e37daf67] [39m[92m+ LibGit2_jll v1.6.4+0[39m
      [90m[29816b5a] [39m[92m+ LibSSH2_jll v1.11.0+1[39m
      [90m[c8ffd9c3] [39m[92m+ MbedTLS_jll v2.28.2+1[39m
      [90m[14a3606d] [39m[92m+ MozillaCACerts_jll v2023.1.10[39m
      [90m[4536629a] [39m[92m+ OpenBLAS_jll v0.3.23+4[39m
      [90m[05823500] [39m[92m+ OpenLibm_jll v0.8.1+2[39m
      [90m[bea87d4a] [39m[92m+ SuiteSparse_jll v7.2.1+1[39m
      [90m[83775a58] [39m[92m+ Zlib_jll v1.2.13+1[39m
      [90m[8e850b90] [39m[92m+ libblastrampoline_jll v5.11.0+0[39m
      [90m[8e850ede] [39m[92m+ nghttp2_jll v1.52.0+1[39m
      [90m[3f19e933] [39m[92m+ p7zip_jll v17.4.0+2[39m
    [36m[1m        Info[22m[39m Packages marked with [32m⌃[39m and [33m⌅[39m have new versions available. Those with [32m⌃[39m may be upgradable, but those with [33m⌅[39m are restricted by compatibility constraints from upgrading. To see why use `status --outdated -m`
    [32m[1m  Activating[22m[39m project at `D:\a\FMIExport.jl\FMIExport.jl\examples`
    

That is all the preperation, that was necessary. Now we can export the FMU. 

TODO The following codecell contains *workardound* code that will need to be modified with the next release


```julia
saveFMU(fmu, fmu_save_path, fmu_source_path; debug=false, compress=false) # feel free to set debug true, disabled for documentation building
#saveFMU(fmu_save_path, fmu_source_path; debug=false, compress=false) this meight be the format after the next release
```

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Generating package ...
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Source package is D:\a\FMIExport.jl\FMIExport.jl\examples\FMI2\BouncingBall, deployed at C:/Users/RUNNER~1/AppData/Local/Temp/fmibuildjl_flOAnu\merged_BouncingBall
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Relative src file path is src\BouncingBall.jl
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... reading FMU template file at C:\Users\runneradmin\.julia\packages\FMIBuild\zfhlh\src/../template/ME/FMU2/src/FMU2_content.jl
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... reading old FMU source file at C:/Users/RUNNER~1/AppData/Local/Temp/fmibuildjl_flOAnu\merged_BouncingBall\src\BouncingBall.jl
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Removing `FMIBUILD_NO_EXPORT_*` blocks ...
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... removed `FMIBUILD_NO_EXPORT_*` blocks.
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Adding/removing dependencies ...
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU]    > Using active environment `D:\a\FMIExport.jl\FMIExport.jl\examples\Project.toml`.
    [32m[1m  Activating[22m[39m project at `D:\a\FMIExport.jl\FMIExport.jl\examples`
    [32m[1m  Activating[22m[39m project at `C:\Users\RUNNER~1\AppData\Local\Temp\fmibuildjl_flOAnu\merged_BouncingBall`
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU]    > Most recent version of `FMIExport` already checked out for FMU, is `D:\a\FMIExport.jl\FMIExport.jl`.
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU]    > Default environment `D:\a\FMIExport.jl\FMIExport.jl\examples\Project.toml` has no dependency on `FMIBase`, adding `FMIBase` from registry.
    

    [32m[1m    Updating[22m[39m registry at `C:\Users\runneradmin\.julia\registries\General.toml`
    [32m[1m   Resolving[22m[39m package versions...
    

    [32m[1m    Updating[22m[39m `C:\Users\runneradmin\AppData\Local\Temp\fmibuildjl_flOAnu\merged_BouncingBall\Project.toml`
      [90m[900ee838] [39m[92m+ FMIBase v1.0.10[39m
    [32m[1m  No Changes[22m[39m to `C:\Users\runneradmin\AppData\Local\Temp\fmibuildjl_flOAnu\merged_BouncingBall\Manifest.toml`
    

    [32m[1m    Updating[22m[39m `C:\Users\runneradmin\AppData\Local\Temp\fmibuildjl_flOAnu\merged_BouncingBall\Project.toml`
      [90m[226f0e26] [39m[91m- FMIBuild v0.3.2[39m
    [32m[1m    Updating[22m[39m `C:\Users\runneradmin\AppData\Local\Temp\fmibuildjl_flOAnu\merged_BouncingBall\Manifest.toml`

    
      [90m[226f0e26] [39m[91m- FMIBuild v0.3.2[39m
      [90m[c27321d9] [39m[91m- Glob v1.3.1[39m
      [90m[9b87118b] [39m[91m- PackageCompiler v2.1.22[39m
      [90m[05181044] [39m[91m- RelocatableFolders v1.0.1[39m
      [90m[6c6a2e73] [39m[91m- Scratch v1.2.1[39m
    [36m[1m        Info[22m[39m We haven't cleaned this depot up for a bit, running Pkg.gc()...
    

    [0mPackageCompiler: bundled libraries:
      ├── Base:
    

    [32m[1m      Active[22m[39m manifest files: 3 found
    [32m[1m      Active[22m[39m artifact files: 9 found
    [32m[1m      Active[22m[39m scratchspaces: 2 found
    [32m[1m     Deleted[22m[39m no artifacts, repos, packages or scratchspaces
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU]    > Removed FMIBuild
    [32m[1m  Activating[22m[39m project at `D:\a\FMIExport.jl\FMIExport.jl\examples`
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... adding/removing dependencies done.
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... generating new FMU source file at C:/Users/RUNNER~1/AppData/Local/Temp/fmibuildjl_flOAnu\merged_BouncingBall\src\BouncingBall.jl
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... generating package done.
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Compiling FMU ...
    

      │    ├── libLLVM-15jl.dll - 116.929 MiB
      │    ├── libatomic-1.dll - 262.531 KiB
      │    ├── libdSFMT.dll - 133.008 KiB
      │    ├── libgcc_s_seh-1.dll - 767.164 KiB
      │    ├── libgfortran-5.dll - 11.162 MiB
      │    ├── libgmp-10.dll - 1.084 MiB
      │    ├── libgmp.dll - 1.085 MiB
      │    ├── libgmpxx-4.dll - 336.805 KiB
      │    ├── libgmpxx.dll - 336.805 KiB
      │    ├── libgomp-1.dll - 1.645 MiB
      │    ├── libjulia-codegen.dll - 103.871 MiB
      │    ├── libjulia-internal.dll - 13.218 MiB
      │    ├── libmpfr-6.dll - 2.518 MiB
      │    ├── libmpfr.dll - 2.519 MiB
      │    ├── libopenlibm.dll - 384.734 KiB
      │    ├── libpcre2-16-0.dll - 711.891 KiB
      │    ├── libpcre2-16.dll - 711.891 KiB
      │    ├── libpcre2-32-0.dll - 683.281 KiB
      │    ├── libpcre2-32.dll - 684.141 KiB
      │    ├── libpcre2-8-0.dll - 774.656 KiB
      │    ├── libpcre2-8.dll - 774.781 KiB
      │    ├── libpcre2-posix-3.dll - 127.047 KiB
      │    ├── libquadmath-0.dll - 1.136 MiB
      │    ├── libssp-0.dll - 143.898 KiB
      │    ├── libstdc++-6.dll - 25.187 MiB
      │    ├── libuv-2.dll - 962.484 KiB
      │    ├── libwinpthread-1.dll - 330.164 KiB
      │    ├── libz.dll - 233.070 KiB
      │    ├── 

    libjulia.dll - 229.312 KiB
      ├── Stdlibs:
      │   ├── LibGit2_jll
      │   │   ├── libgit2.dll - 1.920 MiB
      │   ├── SuiteSparse_jll
      │   │   ├── libamd.dll - 149.055 KiB
      │   │   ├── libbtf.dll - 109.484 KiB
      │   │   ├── libcamd.dll - 154.594 KiB
      │   │   ├── libccolamd.dll - 151.844 KiB
      │   │   ├── libcholmod.dll - 1.666 MiB
      │   │   ├── libcolamd.dll - 134.875 KiB
      │   │   ├── libklu.dll - 322.398 KiB
      │   │   ├── libklu_cholmod.dll - 109.078 KiB
      │   │   ├── libldl.dll - 109.539 KiB
      │   │   ├── librbio.dll - 161.648 KiB
      │   │   ├── libspqr.dll - 621.836 KiB
      │   │   ├── libsuitesparseconfig.dll - 120.320 KiB
      │   │   ├── libumfpack.dll - 1.009 MiB
      │   ├── libblastrampoline_jll
      │   │   ├── libblastrampoline-5.dll - 2.250 MiB
      │   ├── OpenBLAS_jll
      │   │   ├── libopenblas64_.dll - 35.182 MiB
      │   ├── LibSSH2_jll
      │   │   ├── libssh2.dll - 423.625 KiB
      │   ├── LibCURL_jll
      │   │   ├── libcurl-4.dll - 856.172 KiB
      │   ├── MbedTLS_jll
      │   │   ├── libmbedcrypto.dll - 705.938 KiB
      │   │   ├── libmbedtls.dll - 379.828 KiB
      │   │   ├── libmbedx509.dll - 277.492 KiB
      │   ├── nghttp2_jll
      │   │   ├── libnghttp2-14.dll - 868.602 KiB
      Total library file size: 336.291 MiB
    

    [32m[1m Downloading[22m[39m artifact: MKL
    

    [0mPackageCompiler: bundled artifacts:
      ├── IntelOpenMP_jll - 6.198 MiB
      ├── Libiconv_jll - 4.043 MiB
      ├── MKL_jll - 444.082 MiB
    

      ├── OpenSpecFun_jll - 798.680 KiB
      ├── XML2_jll - 11.095 MiB
    

      └── oneTBB_jll - 3.461 MiB
    

      Total artifact file size: 469.659 MiB
    

    [32m[1m Downloading[22m[39m artifact: mingw-w64
    

    - PackageCompiler: creating compiler .ji image (incremental=false)
    

    - PackageCompiler: compiling fresh sysimage (incremental=false)
    

    
    [pid 5768] waiting for IO to finish:
     Handle type        uv_handle_t->data
     timer              0000024e1bb9b5f0->0000024e27f13820


    - PackageCompiler: compiling nonincremental system image

    
    

    [32m[1mPrecompiling[22m[39m

     project...
    

    [32m  ✓ [39m[90mConcreteStructs[39m
    

    [32m  ✓ [39m[90mMuladdMacro[39m
    [32m  ✓ [39m[90mFunctionWrappers[39m
    

    [32m  ✓ [39m[90mExprTools[39m
    

    [32m  ✓ [39m[90mIteratorInterfaceExtensions[39m
    

    [32m  ✓ [39m[90mLazyArtifacts[39m
    

    [32m  ✓ [39m[90mFMICore[39m
    

    [32m  ✓ [39m[90mADTypes[39m
    

    [32m  ✓ [39m[90mOffsetArrays[39m
    

    [32m  ✓ [39m[90mUnPack[39m
    

    [32m  ✓ [39m[90mTest[39m
    

    [32m  ✓ [39m[90mOpenLibm_jll[39m
    

    [32m  ✓ [39m[90mInverseFunctions[39m
    

    [32m  ✓ [39m[90mFillArrays[39m
    

    [32m  ✓ [39m[90mCommonSolve[39m
    

    [32m  ✓ [39m[90mCompilerSupportLibraries_jll[39m
    

    [32m  ✓ [39m[90mSuiteSparse_jll[39m
    

    [32m  ✓ [39m[90mDistributed[39m
    

    [32m  ✓ [39m[90mManualMemory[39m
    

    [32m  ✓ [39m[90mMacroTools[39m
    

    [32m  ✓ [39m[90mPreferences[39m
    

    [32m  ✓ [39m[90mCompat[39m
    

    [32m  ✓ [39m[90mRequires[39m
    

    [32m  ✓ [39m[90mOrderedCollections[39m
    [32m  ✓ [39m[90mDataValueInterfaces[39m
    

    [32m  ✓ [39m[90mReexport[39m
    [32m  ✓ [39m[90mEnumX[39m
    [32m  ✓ [39m[90mSIMDTypes[39m
    

    [32m  ✓ [39m[90mZlib_jll[39m
    

    [32m  ✓ [39m[90mDocStringExtensions[39m
    

    [32m  ✓ [39m[90mIrrationalConstants[39m
    

    [32m  ✓ [39m[90mCompositionsBase[39m
    

    [32m  ✓ [39m[90mCpuId[39m
    

    [32m  ✓ [39m[90mEnzymeCore[39m
    

    [32m  ✓ [39m[90mIfElse[39m
    

    [32m  ✓ [39m[90mConstructionBase[39m
    

    [32m  ✓ [39m[90mCommonWorldInvalidations[39m
    

    [32m  ✓ [39m[90mDataAPI[39m
    

    [32m  ✓ [39m[90mPackageExtensionCompat[39m
    

    [32m  ✓ [39m[90mFastClosures[39m
    

    [32m  ✓ [39m[90mFastLapackInterface[39m
    

    [32m  ✓ [39m[90mFunctors[39m
    

    [32m  ✓ [39m[90mStaticArraysCore[39m

    
    [32m  ✓ [39m[90mFunctionWrappersWrappers[39m
    

    [32m  ✓ [39m[90mRuntimeGeneratedFunctions[39m
    

    [32m  ✓ [39m[90mTableTraits[39m
    

    [32m  ✓ [39m[90mTimerOutputs[39m
    

    [32m  ✓ [39m[90mNaNMath[39m
    

    [32m  ✓ [39m[90mInverseFunctions → InverseFunctionsDatesExt[39m
    

    [32m  ✓ [39m[90mProgressMeter[39m
    

    [32m  ✓ [39m[90mMLStyle[39m
    

    [32m  ✓ [39m[90mThreadingUtilities[39m
    

    [32m  ✓ [39m[90mCommonSubexpressions[39m
    

    [32m  ✓ [39m[90mTruncatedStacktraces[39m
    

    [32m  ✓ [39m[90mJLLWrappers[39m
    

    [32m  ✓ [39m[90mSparseArrays[39m
    [32m  ✓ [39m[90mPrecompileTools[39m
    

    [32m  ✓ [39m[90mCompat → CompatLinearAlgebraExt[39m
    

    [32m  ✓ [39m[90mAdapt[39m
    [32m  ✓ [39m[90mParameters[39m
    

    [32m  ✓ [39m[90mZipFile[39m
    

    [32m  ✓ [39m[90mADTypes → ADTypesEnzymeCoreExt[39m
    

    [32m  ✓ [39m[90mLogExpFunctions[39m
    

    [32m  ✓ [39m[90mConstructionBase → ConstructionBaseLinearAlgebraExt[39m
    

    [32m  ✓ [39m[90mDiffResults[39m
    

    [32m  ✓ [39m[90mTables[39m
    [32m  ✓ [39m[90mInverseFunctions → InverseFunctionsTestExt[39m
    

    [32m  ✓ [39m[90mCompositionsBase → CompositionsBaseInverseFunctionsExt[39m
    

    [32m  ✓ [39m[90mIntelOpenMP_jll[39m
    

    [32m  ✓ [39m[90moneTBB_jll[39m
    

    [32m  ✓ [39m[90mLibiconv_jll[39m
    

    [32m  ✓ [39m[90mOpenSpecFun_jll[39m
    

    [32m  ✓ [39m[90mKLU[39m
    

    [32m  ✓ [39m[90mStatistics[39m
    

    [32m  ✓ [39m[90mFillArrays → FillArraysSparseArraysExt[39m
    

    [32m  ✓ [39m[90mExpronicon[39m
    

    [32m  ✓ [39m[90mRecipesBase[39m
    

    [32m  ✓ [39m[90mStatic[39m
    

    [32m  ✓ [39m[90mChainRulesCore[39m
    

    [32m  ✓ [39m[90mKrylov[39m
    

    [32m  ✓ [39m[90mArrayInterface[39m
    [32m  ✓ [39m[90mDataStructures[39m
    

    [32m  ✓ [39m[90mGPUArraysCore[39m
    

    [32m  ✓ [39m[90mEnzymeCore → AdaptExt[39m
    

    [32m  ✓ [39m[90mOffsetArrays → OffsetArraysAdaptExt[39m
    [32m  ✓ [39m[90mLogExpFunctions → LogExpFunctionsInverseFunctionsExt[39m
    

    [32m  ✓ [39m[90mSetfield[39m
    

    [32m  ✓ [39m[90mMKL_jll[39m
    

    [32m  ✓ [39m[90mAccessors[39m
    

    [32m  ✓ [39m[90mBitTwiddlingConvenienceFunctions[39m
    [32m  ✓ [39m[90mFillArrays → FillArraysStatisticsExt[39m
    

    [32m  ✓ [39m[90mXML2_jll[39m
    

    [32m  ✓ [39m[90mChainRulesCore → ChainRulesCoreSparseArraysExt[39m
    

    [32m  ✓ [39m[90mADTypes → ADTypesChainRulesCoreExt[39m
    

    [32m  ✓ [39m[90mCPUSummary[39m
    

    [32m  ✓ [39m[90mArrayInterface → ArrayInterfaceStaticArraysCoreExt[39m
    [32m  ✓ [39m[90mArrayInterface → ArrayInterfaceSparseArraysExt[39m
    

    [32m  ✓ [39m[90mLogExpFunctions → LogExpFunctionsChainRulesCoreExt[39m
    

    [32m  ✓ [39m[90mArrayInterface → ArrayInterfaceGPUArraysCoreExt[39m
    

    [32m  ✓ [39m[90mSparspak[39m
    [32m  ✓ [39m[90mAccessors → AccessorsDatesExt[39m
    

    [32m  ✓ [39m[90mHostCPUFeatures[39m
    

    [32m  ✓ [39m[90mDifferentiationInterface[39m
    

    [32m  ✓ [39m[90mEzXML[39m
    

    [32m  ✓ [39m[90mSparseMatrixColorings[39m
    

    [32m  ✓ [39m[90mPolyesterWeave[39m
    

    [32m  ✓ [39m[90mSparseConnectivityTracer[39m
    

    [32m  ✓ [39m[90mSciMLStructures[39m
    

    [32m  ✓ [39m[90mSpecialFunctions[39m
    

    [32m  ✓ [39m[90mMaybeInplace[39m
    

    [32m  ✓ [39m[90mFiniteDiff[39m
    

    [32m  ✓ [39m[90mStaticArrayInterface[39m
    

    [32m  ✓ [39m[90mAccessors → AccessorsTestExt[39m
    

    [32m  ✓ [39m[90mDifferentiationInterface → DifferentiationInterfaceChainRulesCoreExt[39m
    

    [32m  ✓ [39m[90mDifferentiationInterface → DifferentiationInterfaceSparseArraysExt[39m
    

    [32m  ✓ [39m[90mSparseConnectivityTracer → SparseConnectivityTracerNaNMathExt[39m
    

    [32m  ✓ [39m[90mSparseConnectivityTracer → SparseConnectivityTracerLogExpFunctionsExt[39m
    

    [32m  ✓ [39m[90mMaybeInplace → MaybeInplaceSparseArraysExt[39m
    

    [32m  ✓ [39m[90mSpecialFunctions → SpecialFunctionsChainRulesCoreExt[39m
    

    [32m  ✓ [39m[90mFiniteDiff → FiniteDiffSparseArraysExt[39m
    

    [32m  ✓ [39m[90mDifferentiationInterface → DifferentiationInterfaceFiniteDiffExt[39m
    

    [32m  ✓ [39m[90mStaticArrayInterface → StaticArrayInterfaceOffsetArraysExt[39m
    

    [32m  ✓ [39m[90mArrayLayouts[39m
    

    [32m  ✓ [39m[90mDifferentiationInterface → DifferentiationInterfaceSparseMatrixColoringsExt[39m
    

    [32m  ✓ [39m[90mSymbolicIndexingInterface[39m
    

    [32m  ✓ [39m[90mSciMLOperators[39m
    

    [32m  ✓ [39m[90mDiffRules[39m
    

    [32m  ✓ [39m[90mCloseOpenIntervals[39m
    

    [32m  ✓ [39m[90mSparseConnectivityTracer → SparseConnectivityTracerSpecialFunctionsExt[39m
    

    [32m  ✓ [39m[90mLayoutPointers[39m
    

    [32m  ✓ [39m[90mArrayLayouts → ArrayLayoutsSparseArraysExt[39m
    

    [32m  ✓ [39m[90mSciMLOperators → SciMLOperatorsStaticArraysCoreExt[39m
    

    [32m  ✓ [39m[90mSciMLOperators → SciMLOperatorsSparseArraysExt[39m
    

    [32m  ✓ [39m[90mRecursiveArrayTools[39m
    

    [32m  ✓ [39m[90mStrideArraysCore[39m
    

    [32m  ✓ [39m[90mForwardDiff[39m
    

    [32m  ✓ [39m[90mLazyArrays[39m
    [32m  ✓ [39m[90mRecursiveArrayTools → RecursiveArrayToolsSparseArraysExt[39m
    

    [32m  ✓ [39m[90mPolyester[39m
    

    [32m  ✓ [39m[90mPreallocationTools[39m
    

    [32m  ✓ [39m[90mNLSolversBase[39m
    

    [32m  ✓ [39m[90mVectorizationBase[39m
    

    [32m  ✓ [39m[90mDifferentiationInterface → DifferentiationInterfaceForwardDiffExt[39m
    

    [32m  ✓ [39m[90mFastBroadcast[39m
    

    [32m  ✓ [39m[90mRecursiveArrayTools → RecursiveArrayToolsForwardDiffExt[39m
    

    [32m  ✓ [39m[90mLineSearches[39m
    

    [32m  ✓ [39m[90mSLEEFPirates[39m
    

    [32m  ✓ [39m[90mRecursiveArrayTools → RecursiveArrayToolsFastBroadcastExt[39m
    

    [32m  ✓ [39m[90mSciMLBase[39m
    

    [32m  ✓ [39m[90mSciMLBase → SciMLBaseChainRulesCoreExt[39m
    

    [32m  ✓ [39m[90mLoopVectorization[39m
    

    [32m  ✓ [39m[90mLoopVectorization → SpecialFunctionsExt[39m
    

    [32m  ✓ [39m[90mSciMLJacobianOperators[39m
    

    [32m  ✓ [39m[90mLoopVectorization → ForwardDiffExt[39m
    

    [32m  ✓ [39m[90mDiffEqBase[39m
    

    [32m  ✓ [39m[90mLineSearch[39m
    

    [32m  ✓ [39m[90mDiffEqBase → DiffEqBaseSparseArraysExt[39m
    

    [32m  ✓ [39m[90mTriangularSolve[39m
    

    [32m  ✓ [39m[90mLineSearch → LineSearchLineSearchesExt[39m
    

    [32m  ✓ [39m[90mDiffEqBase → DiffEqBaseChainRulesCoreExt[39m
    

    [32m  ✓ [39m[90mRecursiveFactorization[39m
    

    [32m  ✓ [39m

    [90mSimpleNonlinearSolve[39m
    

    [32m  ✓ [39m[90mSimpleNonlinearSolve → SimpleNonlinearSolveChainRulesCoreExt[39m
    

    [32m  ✓ [39m[90mLinearSolve[39m
    

    [32m  ✓ [39m[90mLinearSolve → LinearSolveRecursiveArrayToolsExt[39m
    

    [32m  ✓ [39m[90mLinearSolve → LinearSolveEnzymeExt[39m
    

    [32m  ✓ [39m[90mNonlinearSolve[39m
    

    [32m  ✓ [39m[90mDiffEqCallbacks[39m
    

    [32m  ✓ [39mFMIBase
    

    [32m  ✓ [39m[90mFMIBase → ForwardDiffExt[39m
    

    [32m  ✓ [39mFMIExport
    

    [32m  ✓ [39mBouncingBall
    

      172 dependencies successfully precompiled in 174 seconds
    

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39mPackageCompiler: Executing C:\Users\RUNNER~1\AppData\Local\Temp\fmibuildjl_flOAnu\merged_BouncingBall\src\BouncingBall.jl => C:\Users\RUNNER~1\AppData\Local\Temp\jl_packagecompiler_DZ2HDN\jl_93C6.tmp
    

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39mPackageCompiler: Done
    - PackageCompiler: compiling incremental system image
    

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... compiling FMU done.
    

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Patching libjulia.dll @ `C:\Users\RUNNER~1\AppData\Local\Temp\fmibuildjl_flOAnu\BouncingBall\binaries\win64`...
    

    [36m[1m┌ [22m[39m[36m[1mInfo: [22m[39mFound embedded libpath
    [36m[1m│ [22m[39m  libpath =
    [36m[1m│ [22m[39m   6-element Vector{SubString{String}}:
    [36m[1m│ [22m[39m    "../bin/libgcc_s_seh-1.dll"
    [36m[1m│ [22m[39m    "../bin/libopenlibm.dll"
    [36m[1m│ [22m[39m    "@"
    [36m[1m│ [22m[39m    "@../bin/libjulia-internal.dll"
    [36m[1m│ [22m[39m    "@../bin/libjulia-codegen.dll"
    [36m[1m│ [22m[39m    ""
    [36m[1m└ [22m[39m  libpath_offset = 14849
    

    	modelDescription.xml
    

    [36m[1m┌ [22m[39m[36m[1mInfo: [22m[39mFiltered libpath
    [36m[1m│ [22m[39m  libpath =
    [36m[1m│ [22m[39m   6-element Vector{AbstractString}:
    [36m[1m│ [22m[39m    "libgcc_s_seh-1.dll"
    [36m[1m│ [22m[39m    "libopenlibm.dll"
    [36m[1m│ [22m[39m    "@"
    [36m[1m│ [22m[39m    "@libjulia-internal.dll"
    [36m[1m│ [22m[39m    "@libjulia-codegen.dll"
    [36m[1m└ [22m[39m    ""
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... patching libjulia.dll done.
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Building model description ...
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... building model description done.
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Zipping FMU ...
    

    	binaries/include/FMU2_init.h
    	binaries/include/julia_init.h
    	binaries/share/julia/LocalPreferences.toml
    	binaries/share/julia/Project.toml
    	binaries/share/julia/cert.pem
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/libimalloc.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_avx2.2.dll
    

    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_avx512.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_blacs_ilp64.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_blacs_intelmpi_ilp64.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_blacs_intelmpi_lp64.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_blacs_lp64.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_blacs_msmpi_ilp64.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_blacs_msmpi_lp64.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_cdft_core.2.dll

    
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_core.2.dll
    

    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_def.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_intel_thread.2.dll
    

    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_mc3.2.dll
    

    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_pgi_thread.2.dll
    

    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_rt.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_scalapack_ilp64.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_scalapack_lp64.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_sequential.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_tbb_thread.2.dll
    

    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_vml_avx2.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_vml_avx512.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_vml_cmpt.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_vml_def.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/bin/mkl_vml_mc3.2.dll
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/share/licenses/MKL/license.txt
    	binaries/share/julia/artifacts/07d4e2a4c926ce5a99f7ff402617dc9caa2187a0/share/licenses/MKL/tpp.txt
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/bin/libiomp5md.dll
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/bin/libiomp5md_db.dll
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/bin/libiompstubs5md.dll
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/bin/omptarget.dll
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/bin/omptarget.rtl.level0.dll
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/bin/omptarget.rtl.opencl.dll
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/bin/omptarget.rtl.unified_runtime.dll
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/bin/omptarget.sycl.wrap.dll
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/lib/libiomp5md.lib
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/lib/libiompstubs5md.lib
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/share/licenses/IntelOpenMP/LICENSE.txt
    	binaries/share/julia/artifacts/3946a6b7e5b0c7955b4a4e68280c32b229774d34/share/licenses/IntelOpenMP/third-party-programs.txt
    	binaries/share/julia/artifacts/3e683ec5ca945a5aca74c49e8cccdf37c19b84a3/bin/libopenspecfun.dll
    	binaries/share/julia/artifacts/3e683ec5ca945a5aca74c49e8cccdf37c19b84a3/include/Faddeeva.h
    	binaries/share/julia/artifacts/3e683ec5ca945a5aca74c49e8cccdf37c19b84a3/lib/libopenspecfun.a
    	binaries/share/julia/artifacts/3e683ec5ca945a5aca74c49e8cccdf37c19b84a3/logs/OpenSpecFun/OpenSpecFun.log.gz
    	binaries/share/julia/artifacts/3e683ec5ca945a5aca74c49e8cccdf37c19b84a3/share/licenses/OpenSpecFun/LICENSE.md
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/bin/libxml2-2.dll
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/bin/xml2-config
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/bin/xmlcatalog.exe
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/bin/xmllint.exe
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/HTMLparser.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/HTMLtree.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/SAX.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/SAX2.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/c14n.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/catalog.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/chvalid.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/debugXML.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/dict.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/encoding.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/entities.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/globals.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/hash.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/list.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/nanoftp.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/nanohttp.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/parser.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/parserInternals.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/pattern.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/relaxng.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/schemasInternals.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/schematron.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/threads.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/tree.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/uri.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/valid.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xinclude.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xlink.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlIO.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlautomata.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlerror.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlexports.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlmemory.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlmodule.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlreader.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlregexp.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlsave.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlschemas.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlschemastypes.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlstring.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlunicode.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlversion.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xmlwriter.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xpath.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xpathInternals.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/include/libxml2/libxml/xpointer.h
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/lib/libxml2.dll.a
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/lib/cmake/libxml2/libxml2-config.cmake
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/lib/pkgconfig/libxml-2.0.pc
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/aclocal/libxml.m4
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/general.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/home.png
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/index.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/left.png
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-HTMLparser.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-HTMLtree.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-SAX.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-SAX2.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-c14n.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-catalog.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-chvalid.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-debugXML.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-dict.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-encoding.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-entities.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-globals.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-hash.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-list.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-nanoftp.html

    
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-nanohttp.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-parser.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-parserInternals.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-pattern.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-relaxng.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-schemasInternals.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-schematron.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-threads.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-tree.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-uri.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-valid.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xinclude.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xlink.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlIO.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlautomata.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlerror.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlexports.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlmemory.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlmodule.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlreader.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlregexp.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlsave.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlschemas.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlschemastypes.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlstring.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlunicode.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlversion.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xmlwriter.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xpath.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xpathInternals.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2-xpointer.html
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/libxml2.devhelp2
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/right.png
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/style.css
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/gtk-doc/html/libxml2/up.png
    	binaries/share/julia/artifacts/55b1a3d509033b500e2c33e984a35f8ed481879a/share/licenses/XML2/Copyright
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/bin/iconv.exe
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/bin/libcharset-1.dll

    
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/bin/libiconv-2.dll
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/include/iconv.h
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/include/libcharset.h
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/include/localcharset.h
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/lib/libcharset.a
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/lib/libcharset.dll.a
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/lib/libiconv.a
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/lib/libiconv.dll.a
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/lib/pkgconfig/iconv.pc
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/doc/iconv.1.html
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/doc/iconv.3.html
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/doc/iconv_close.3.html
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/doc/iconv_open.3.html
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/doc/iconv_open_into.3.html
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/doc/iconvctl.3.html
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/licenses/Libiconv/COPYING
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/man/man1/iconv.1
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/man/man3/iconv.3
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/man/man3/iconv_close.3
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/man/man3/iconv_open.3
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/man/man3/iconv_open_into.3
    	binaries/share/julia/artifacts/79e4bc6534ea5a11e42eaf15947a2272949e4865/share/man/man3/iconvctl.3
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/bin/libtbb12.dll
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/bin/libtbbmalloc.dll
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/bin/libtbbmalloc_proxy.dll
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/blocked_range.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/blocked_range2d.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/blocked_range3d.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/blocked_rangeNd.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/cache_aligned_allocator.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/collaborative_call_once.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/combinable.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_hash_map.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_lru_cache.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_map.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_priority_queue.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_queue.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_set.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_unordered_map.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_unordered_set.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/concurrent_vector.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/enumerable_thread_specific.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/flow_graph.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/flow_graph_abstractions.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/global_control.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/info.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/memory_pool.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/null_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/null_rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/parallel_for.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/parallel_for_each.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/parallel_invoke.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/parallel_pipeline.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/parallel_reduce.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/parallel_scan.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/parallel_sort.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/partitioner.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/profiling.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/queuing_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/queuing_rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/scalable_allocator.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/spin_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/spin_rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/task.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/task_arena.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/task_group.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/task_scheduler_observer.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/tbb_allocator.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/tbbmalloc_proxy.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/tick_count.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/version.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_aggregator.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_aligned_space.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_allocator_traits.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_assert.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_attach.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_concurrent_queue_base.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_concurrent_skip_list.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_concurrent_unordered_base.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_config.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_containers_helpers.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_exception.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_export.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_body_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_cache_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_indexer_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_item_buffer_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_join_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_node_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_node_set_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_nodes_deduction.h

    
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_tagged_buffer_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_trace_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_flow_graph_types_impl.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_hash_compare.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_intrusive_list_node.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_machine.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_mutex_common.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_namespace_injection.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_node_handle.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_pipeline_filters.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_pipeline_filters_deduction.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_range_common.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_rtm_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_rtm_rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_scoped_lock.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_segment_table.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_small_object_pool.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_string_resource.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_task.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_task_handle.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_template_helpers.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_utils.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/oneapi/tbb/detail/_waitable_atomic.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/blocked_range.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/blocked_range2d.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/blocked_range3d.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/blocked_rangeNd.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/cache_aligned_allocator.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/collaborative_call_once.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/combinable.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_hash_map.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_lru_cache.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_map.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_priority_queue.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_queue.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_set.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_unordered_map.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_unordered_set.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/concurrent_vector.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/enumerable_thread_specific.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/flow_graph.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/flow_graph_abstractions.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/global_control.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/info.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/memory_pool.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/null_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/null_rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/parallel_for.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/parallel_for_each.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/parallel_invoke.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/parallel_pipeline.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/parallel_reduce.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/parallel_scan.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/parallel_sort.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/partitioner.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/profiling.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/queuing_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/queuing_rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/scalable_allocator.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/spin_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/spin_rw_mutex.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/task.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/task_arena.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/task_group.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/task_scheduler_observer.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/tbb.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/tbb_allocator.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/tbbmalloc_proxy.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/tick_count.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/include/tbb/version.h
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/lib/libtbb12.dll.a
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/lib/libtbbmalloc.dll.a
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/lib/libtbbmalloc_proxy.dll.a
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/lib/cmake/TBB/TBBConfig.cmake
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/lib/cmake/TBB/TBBConfigVersion.cmake
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/lib/cmake/TBB/TBBTargets-release.cmake
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/lib/cmake/TBB/TBBTargets.cmake
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/lib/pkgconfig/tbb.pc
    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/share/doc/TBB/README.md
    

    	binaries/share/julia/artifacts/a80da89d0fb528913efffff81a0cf1e2ab988cb2/share/licenses/oneTBB/LICENSE.txt
    	binaries/win64/BouncingBall.dll

    
    	binaries/win64/libLLVM-15jl.dll
    

    	binaries/win64/libamd.dll
    	binaries/win64/libatomic-1.dll
    	binaries/win64/libblastrampoline-5.dll
    	binaries/win64/libbtf.dll
    	binaries/win64/libcamd.dll
    	binaries/win64/libccolamd.dll
    	binaries/win64/libcholmod.dll
    	binaries/win64/libcolamd.dll
    	binaries/win64/libcurl-4.dll
    	binaries/win64/libdSFMT.dll
    	binaries/win64/libgcc_s_seh-1.dll
    	binaries/win64/libgfortran-5.dll
    	binaries/win64/libgit2.dll
    	binaries/win64/libgmp-10.dll
    	binaries/win64/libgmp.dll
    	binaries/win64/libgmpxx-4.dll
    	binaries/win64/libgmpxx.dll
    	binaries/win64/libgomp-1.dll
    	binaries/win64/libjulia-codegen.dll

    
    	binaries/win64/libjulia-internal.dll
    	binaries/win64/libjulia.dll
    	binaries/win64/libklu.dll
    	binaries/win64/libklu_cholmod.dll
    	binaries/win64/libldl.dll
    	binaries/win64/libmbedcrypto.dll
    	binaries/win64/libmbedtls.dll
    	binaries/win64/libmbedx509.dll
    	binaries/win64/libmpfr-6.dll
    	binaries/win64/libmpfr.dll
    	binaries/win64/libnghttp2-14.dll
    	binaries/win64/libopenblas64_.dll
    	binaries/win64/libopenlibm.dll
    

    	binaries/win64/libpcre2-16-0.dll
    	binaries/win64/libpcre2-16.dll
    	binaries/win64/libpcre2-32-0.dll
    	binaries/win64/libpcre2-32.dll
    	binaries/win64/libpcre2-8-0.dll
    	binaries/win64/libpcre2-8.dll
    	binaries/win64/libpcre2-posix-3.dll
    	binaries/win64/libquadmath-0.dll
    	binaries/win64/librbio.dll
    	binaries/win64/libspqr.dll
    	binaries/win64/libssh2.dll
    	binaries/win64/libssp-0.dll
    	binaries/win64/libstdc++-6.dll
    

    	binaries/win64/libsuitesparseconfig.dll
    	binaries/win64/libumfpack.dll
    	binaries/win64/libuv-2.dll
    	binaries/win64/libwinpthread-1.dll
    	binaries/win64/libz.dll
    

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... zipping FMU done.
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] Clean up ...
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39m[Build FMU] ... clean up done.
    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39mFMU-Export succeeded after 15m 59s (1.0% packing time)
    




    true



Now we will grab the generated FMU and move it to a path, where it will be included in this documentation


```julia
mkpath("Export-BouncingBall_files")
cp(fmu_save_path, joinpath("Export-BouncingBall_files", "BouncingBall.fmu"))
```




    "Export-BouncingBall_files\\BouncingBall.fmu"



One current limitation of Julia-FMUs is, that they can not be imported back into Julia, as it is currently not allowed having two Julia-sys-images existing at the same time within the same process. (The Julia FMU comes bundeled with its own image) 

TODO Therefore we will test our generated FMU in Python unsing FMPy.
