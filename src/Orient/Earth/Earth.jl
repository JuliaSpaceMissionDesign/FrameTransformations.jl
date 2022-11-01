include("types.jl")
include("fundamentals.jl")

function build_cio_series(fname::Symbol, iau_model::Symbol, 
                          cpoly::AbstractVector{<:Number}, 
                          ctrig::AbstractVector{<:AbstractVector{<:IAUSeries}})

    # Generates ad-hoc functions to compute the IAU Series including only 
    # the non-null coefficients. 

    # Maps IAUSeries coefficients to Fundamental Argument name
    MAPPING = SVector(:Mₐ, :Sₐ, :uₘ, :Dₛ, :Ωₘ, :λ_Me, :λ_Ve, :λ_Ea, 
                      :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ)

    # cpoly and ctrig must have the same number of elements
    nblocks = length(cpoly)

    fcn_body = Expr(:block)

    # Final @evalpoly call expression (depends on nblocks)
    poly_call = Expr(:macrocall, Symbol("@evalpoly"), :, :t)
    
    # Stores the indexes of all the contributions that have already been computed 
    dct = Dict{Int, Vector{Int}}()

    # Assings all the polynomial contributions
    for i = 1:nblocks 
        dct[i] = Vector{Int}() # initialises dictionary keys 
        push!(fcn_body.args, Expr(:(=), Symbol("w$i"), cpoly[i]))
        push!(poly_call.args, Symbol("w$i"))
    end

    for i = 1:nblocks 

        # Parses the trigonometric series
        for s in ctrig[i]
            
            # Computes the ARGUMENT expression
            argv = Vector{Expr}()
            for (j, n) in enumerate(s.N) 
                if n != 0 
                    f = Expr(:(.), :fa, QuoteNode(MAPPING[j]))
                    push!(argv, Expr(:call, :(*), n, f))
                end
            end

            arge = length(argv) > 0 ? Expr(:call, :(+), argv...) : 0.

            # Stores current block expression
            blkₑ = Vector{Expr}()

            # Re-parses the entire trigonometric coefficients vector to check 
            # for elements with same arg coeffs 
            cn, sn = false, false 
            
            for ii = 1:nblocks 
                fcc, scc = 0., 0. 

                for (k, sk) in enumerate(ctrig[ii])
                    if sk.N == s.N && !(k in dct[ii])
                        # Stores the current index
                        push!(dct[ii], k)
                        
                        # Groups together all the sin\cos terms
                        fcc += sk.fc
                        scc += sk.sc 
   
                        # Stores whether sin\cos are actually computed
                        sn = sn ? sn : sk.fc != 0.
                        cn = cn ? cn : sk.sc != 0.
                    end
                end

                sumₑ = Expr(:call, :(+))
                fcc != 0. && push!(sumₑ.args, Expr(:call, :(*), fcc, :sarg))
                scc != 0. && push!(sumₑ.args, Expr(:call, :(*), scc, :carg))

                if fcc != 0. || scc != 0. 
                    push!(blkₑ, Expr(:(+=), Symbol("w$ii"), sumₑ))    
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

    push!(fcn_body.args, poly_call)

    # Creates function call Expression
    fcall = Expr(:call, fname)
    push!(fcall.args, Expr(:(::), iau_model)) # Adds model
    push!(fcall.args, Expr(:(::), :t, :Number)) # Adds time 
    push!(fcall.args, Expr(:(::), :fa, :FundamentalArguments)) # Adds FA

    # Assembles Function head and body
    fcn = Expr(:function, fcall, fcn_body)
    return eval(fcn)
end

function build_nutation_series(fname::Symbol, iau_model::Symbol, 
                               ψseries::AbstractVector{N}, 
                               ϵseries::AbstractVector{N}) where {N <: AbstractVector{<:IAUSeries} }

    # Generates ad-hoc functions to compute the IAU Series including only 
    # the non-null coefficients. 

    # Maps IAUSeries coefficients to Fundamental Argument name
    MAPPING = SVector(:Mₐ, :Sₐ, :uₘ, :Dₛ, :Ωₘ, :λ_Me, :λ_Ve, :λ_Ea, 
                      :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ)

    # ψseries and ϵseries must have the same number of elements (should be 2)
    nblocks = length(ψseries)

    fcn_body = Expr(:block)

    # Final @evalpoly call expression (depends on nblocks)
    pcall = Vector{Expr}(undef, 2); 
    for i = 1:2 
        var = i == 1 ? :dψ : :dϵ

        # Automatically includes transformation to radians
        pcall[i] = Expr(:(=), var, Expr(:call, :(/), Expr(:call, :(*), 
                            Expr(:macrocall, Symbol("@evalpoly"), :, :t), 
                        π*1e-7), 648000))
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
                argv = Vector{Expr}()
                for (j, n) in enumerate(s.N) 
                    if n != 0 
                        f = Expr(:(.), :fa, QuoteNode(MAPPING[j]))
                        push!(argv, Expr(:call, :(*), n, f))
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
                        fcc, scc = 0., 0. 

                        # Parses the inner series 
                        for (k, sk) in enumerate(vct[ii])
                            if sk.N == s.N && !(k in dct[(v, ii)])

                                # Stores the current index
                                push!(dct[(v, ii)], k)
                                
                                # Groups together all the sin\cos terms
                                fcc += sk.fc 
                                scc += sk.sc
        
                                # Stores whether sin\cos are actually computed
                                # For ψ fc = sin, sc = cos, for ϵ fc = cos, sc = sin
                                cn = cn ? cn : (v == 1 ? sk.sc : sk.fc) != 0.
                                sn = sn ? sn : (v == 1 ? sk.fc : sk.sc) != 0.
                            end
                        end

                        sumₑ = Expr(:call, :(+))

                        if fcc != 0. 
                            push!(sumₑ.args, 
                                Expr(:call, :(*), fcc, v == 1 ? (:sarg) : (:carg)))
                        end

                        if scc != 0. 
                            push!(sumₑ.args, 
                                Expr(:call, :(*), scc, v == 1 ? (:carg) : (:sarg)))
                        end

                        if fcc != 0. || scc != 0. 
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
    
    res = Expr(:call, :(./), Expr(:call, :(.*), Expr(:tuple, :dψ, :dϵ), 
                                 π*1e-7), 648000)

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

include("iers.jl")
include("nutation.jl")

