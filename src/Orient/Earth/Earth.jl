include("types.jl")
include("fundamentals.jl")
include("eop.jl")

function build_cio_series(fname::Symbol, iau_model::Symbol, 
    cpoly::AbstractVector{<:AbstractVector{<:Number}}, 
    ctrig::AbstractVector{<:AbstractVector{<:AbstractVector{<:IAUSeries}}})

    nseries = length(cpoly)
    if nseries != length(ctrig) 
        throw(ArgumentError("[Orient] cpoly and ctrig must have same number of series."))
    end 

    # Generates ad-hoc functions to compute the IAU Series including only 
    # the non-null coefficients. 

    # The elements of cpoly must always be >= ctrig
    fcn_body = Expr(:block)

    # Maps IAUSeries coefficients to Fundamental Argument name
    MAPPING = SVector(:Mₐ, :Sₐ, :uₘ, :Dₛ, :Ωₘ, :λ_Me, :λ_Ve, :λ_Ea, 
                      :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ)
    
    # Pre-computes all the Fundamental Arguments multiplications
    fad = Dict{Symbol, Vector{Int}}()
    for (m, fm) in enumerate(MAPPING) 
        fad[fm] = Vector{Int}()

        # Finds all the multiplicative factors of argument `fm`
        for i in eachindex(ctrig)
            for j in eachindex(ctrig[i])
                for s in ctrig[i][j]
                    if !(s.N[m] in fad[fm]) && s.N[m] != 0 
                        push!(fad[fm], s.N[m])
                    end
                end
            end
        end

        # Register the expression for that facotr
        for (k, mk) in enumerate(fad[fm])
            fae = Expr(:(.), :fa, QuoteNode(fm))
            push!(fcn_body.args, Expr(:(=), Symbol("$fm$k"), Expr(:call, :(*), mk, fae)))
        end
    end

    # Final @evalpoly call expression (depends on nblocks)
    pcall = Vector{Expr}(undef, nseries)
    for j = 1:nseries 
        var = Symbol("x$(j)x")

        # Automatically includes transformation to radians
        pcall[j] = Expr(:(=), var, Expr(:call, :(/), Expr(:call, :(*), 
                        Expr(:macrocall, Symbol("@evalpoly"), :, :t), π*1e-6), 648000))

    end

    # Stores the indexes of all the contributions that have already been computed 
    dct = Dict{Tuple{Int, Int}, Vector{Int}}()

    for j = 1:nseries 
        # Assings all the polynomial contributions
        for i in eachindex(cpoly[j])
            dct[(j, i)] = Vector{Int}() # initialises dictionary keys 
            push!(fcn_body.args, Expr(:(=), Symbol("x$(j)x$i"), cpoly[j][i]))
            push!(pcall[j].args[2].args[2].args[2].args, Symbol("x$(j)x$i"))
        end
    end

    for jj = 1:nseries 
        for i in eachindex(ctrig[jj]) 
            # Parses the trigonometric series
            for s in ctrig[jj][i]

                # Computes the ARGUMENT expression
                argv = Vector{Symbol}()
                for (j, n) in enumerate(s.N) 
                    if n != 0 
                        push!(argv, Symbol("$(MAPPING[j])$(argmin(abs.(fad[MAPPING[j]] .- n)))"))
                    end
                end

                arge = length(argv) > 0 ? Expr(:call, :(+), argv...) : 0.

                # Stores current block expression
                blkₑ = Vector{Expr}()

                # Re-parses the entire trigonometric coefficients vector to check 
                # for elements with same arg coeffs 
                cn, sn = false, false 

                for jjj = 1:nseries 
                    for ii in eachindex(ctrig[jjj]) 
                        sc, cc = 0., 0. 

                        for (k, sk) in enumerate(ctrig[jjj][ii])
                            if sk.N == s.N && !(k in dct[(jjj, ii)])
                                # Stores the current index
                                push!(dct[(jjj, ii)], k)
                                
                                # Groups together all the sin\cos terms
                                sc += sk.sc
                                cc += sk.cc 

                                # Stores whether sin\cos are actually computed
                                sn = sn ? sn : sk.sc != 0.
                                cn = cn ? cn : sk.cc != 0.
                            end
                        end

                        sumₑ = Expr(:call, :(+))
                        sc != 0. && push!(sumₑ.args, Expr(:call, :(*), sc, :sarg))
                        cc != 0. && push!(sumₑ.args, Expr(:call, :(*), cc, :carg))

                        if sc != 0. || cc != 0. 
                            push!(blkₑ, Expr(:(+=), Symbol("x$(jjj)x$ii"), sumₑ))    
                        end      
                    end
                end

                # Adds all the expressions 
                if cn || sn 

                    # Reassigns the arg variable
                    push!(fcn_body.args, Expr(:(=), :arg, arge))

                    # Computes and stores sin(ARG) and cos(ARG)
                    cn && push!(fcn_body.args, Expr(:(=), :carg, Expr(:call, :cos, :arg)))
                    sn && push!(fcn_body.args, Expr(:(=), :sarg, Expr(:call, :sin, :arg)))
                        
                    push!(fcn_body.args, blkₑ...)
                end
            end
        end
    end

    for j = 1:nseries 
        push!(fcn_body.args, pcall[j])
    end

    # Adds return types
    if nseries > 1 
        rtupl= Expr(:tuple)
        for j = 1:nseries 
            push!(rtupl.args, Symbol("x$(j)x"))
        end

        push!(fcn_body.args, Expr(:return, rtupl))
    end

    # Creates function call Expression
    fcall = Expr(:call, fname)
    push!(fcall.args, Expr(:(::), iau_model)) # Adds model
    push!(fcall.args, Expr(:(::), :t, :Number)) # Adds time 
    push!(fcall.args, Expr(:(::), :fa, :FundamentalArguments)) # Adds FA

    # Assembles Function head and body
    fcn = Expr(:function, fcall, fcn_body)
    return eval(fcn)
end

function build_nutation_series(fname::Symbol, iau_model::Symbol, ψseries::AbstractVector{N}, 
            ϵseries::AbstractVector{N}) where {N <: AbstractVector{<:IAUSeries} }

    # Generates ad-hoc functions to compute the IAU Series including only 
    # the non-null coefficients. 

    # Maps IAUSeries coefficients to Fundamental Argument name
    MAPPING = SVector(:Mₐ, :Sₐ, :uₘ, :Dₛ, :Ωₘ, :λ_Me, :λ_Ve, :λ_Ea, 
                      :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ)

    # ψseries and ϵseries must have the same number of elements (should be 2)
    nblocks = length(ψseries)
    fcn_body = Expr(:block)

    # Pre-computes all the Fundamental Arguments multiplications
    fad = Dict{Symbol, Vector{Int}}()
    for (m, fm) in enumerate(MAPPING) 
        fad[fm] = Vector{Int}()

        # Finds all the multiplicative factors of argument `fm`
        for i = 1:nblocks 
            for series in (ψseries, ϵseries)
                for s in series[i]
                    if !(s.N[m] in fad[fm]) && s.N[m] != 0 
                        push!(fad[fm], s.N[m])
                    end
                end
            end
        end

        # Register the expression for that facotr
        for (k, mk) in enumerate(fad[fm])
            fae = Expr(:(.), :fa, QuoteNode(fm))
            push!(fcn_body.args, Expr(:(=), Symbol("$fm$k"), Expr(:call, :(*), mk, fae)))
        end
    end

    # Final @evalpoly call expression (depends on nblocks)
    pcall = Vector{Expr}(undef, 2); 
    for i = 1:2 
        var = i == 1 ? :dψ : :dϵ

        # Automatically includes transformation to radians
        pcall[i] = Expr(:(=), var, Expr(:call, :(/), Expr(:call, :(*), 
                        Expr(:macrocall, Symbol("@evalpoly"), :, :t), π*1e-6), 648000))
    end

    # Stores the indexes of all the contributions that have already been computed 
    # First Tuple index refers to ψ or ϵ, the second to the corresponding block
    dct = Dict{Tuple{Int, Int}, Vector{Int}}()

    # Initialises the variables 
    for i = 1:nblocks 
        for j = 1:2 # Adds ψ and ϵ contributions 
            dct[(j, i)] = Vector{Int}() # initialises dictionary keys 

            var = j == 1 ? Symbol("dψ$(i-1)") : Symbol("dϵ$(i-1)") 
            push!(fcn_body.args, Expr(:(=), var, 0.0))
            push!(pcall[j].args[2].args[2].args[2].args, var)
        end
    end

    for i = 1:nblocks 
        for series in (ψseries, ϵseries)

            # Parses ψ or ϵ block 
            for s in series[i]

                # Computes the ARGUMENT expression
                argv = Vector{Symbol}()
                for (j, n) in enumerate(s.N) 
                    if n != 0 
                        push!(argv, Symbol("$(MAPPING[j])$(argmin(abs.(fad[MAPPING[j]] .- n)))"))
                    end
                end

                arge = length(argv) > 0 ? Expr(:call, :(+), argv...) : 0.

                # Stores current block expression
                blkₑ = Vector{Expr}()

                # Re-parses the entire ψ and ϵ vectors to check 
                # for elements with same arg coeffs 
                cn, sn = false, false 

                for ii = 1:nblocks 

                    # Parse ψ and ϵ
                    for (v, vct) in enumerate([ψseries, ϵseries])
                        sc, cc = 0., 0. 

                        # Parses the inner series 
                        for (k, sk) in enumerate(vct[ii])
                            if sk.N == s.N && !(k in dct[(v, ii)])

                                # Stores the current index
                                push!(dct[(v, ii)], k)
                                
                                # Groups together all the sin\cos terms
                                sc += sk.sc 
                                cc += sk.cc

                                # Stores whether sin\cos are actually computed
                                sn = sn ? sn : sk.sc != 0.
                                cn = cn ? cn : sk.cc != 0.
                            end
                        end

                        sumₑ = Expr(:call, :(+))
                        sc != 0. && push!(sumₑ.args, Expr(:call, :(*), sc, :sarg))
                        cc != 0. && push!(sumₑ.args, Expr(:call, :(*), cc, :carg))

                        if sc != 0. || cc != 0. 
                            var = v == 1 ? "ψ" : "ϵ"
                            push!(blkₑ, Expr(:(+=), Symbol("d"*var*"$(ii-1)"), sumₑ))
                        end
                    end
                end

                # Adds all the expressions 
                if cn || sn 

                    # Reassigns the arg variable
                    push!(fcn_body.args, Expr(:(=), :arg, arge))

                    # Computes and stores sin(ARG) and cos(ARG)
                    cn && push!(fcn_body.args, Expr(:(=), :carg, Expr(:call, :cos, :arg)))
                    sn && push!(fcn_body.args, Expr(:(=), :sarg, Expr(:call, :sin, :arg)))
                        
                    push!(fcn_body.args, blkₑ...)
                end
            end
        end
    end

    for i = 1:2 
        push!(fcn_body.args, pcall[i])
    end

    # Adds return types 
    push!(fcn_body.args, Expr(:return, Expr(:tuple, :dψ, :dϵ)))

    # Creates function call Expression
    fcall = Expr(:call, fname)
    push!(fcall.args, Expr(:(::), iau_model)) # Adds model
    push!(fcall.args, Expr(:(::), :t, :Number)) # Adds time 
    push!(fcall.args, Expr(:(::), :fa, :FundamentalArguments)) # Adds FA

    # Assembles Function head and body
    fcn = Expr(:function, fcall, fcn_body)
    return eval(fcn)
end

include("obliquity.jl")
include("nutation.jl")
include("precession.jl")

include("iers.jl")
