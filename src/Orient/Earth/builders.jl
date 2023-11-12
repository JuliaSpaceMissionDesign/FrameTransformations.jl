
# These functions are used to automatically generated optimised functions that 
# compute the series associated to the IAU constants.

# Used to map the IAUSeries coefficients to the fundamental argument name 
const _FA_MAPPING = SVector(
    :Mₐ, :Sₐ, :uₘ, :Dₛ, :Ωₘ, :λ_Me, :λ_Ve, :λ_Ea, :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ
)

function _assemble_function(iau_model, fname, fbody)

    # Creates function call Expression
    fcall = Expr(:call, fname)

    push!(fcall.args, Expr(:(::), iau_model))   # Adds model
    push!(fcall.args, Expr(:(::), :t, :Number)) # Adds time 
    push!(fcall.args, Expr(:(::), :fa, :FundamentalArguments)) # Adds FA

    # Assembles Function head and body
    return Expr(:function, fcall, fbody)

end

function build_cio_series(
    fname::Symbol,
    iau_model::Symbol,
    cpoly::AbstractVector{<:AbstractVector{<:Number}},
    ctrig::AbstractVector{<:AbstractVector{<:AbstractVector{<:IAUSeries}}},
)
    nseries = length(cpoly)
    if nseries != length(ctrig)
        throw(ArgumentError("[Orient] cpoly and ctrig must have same number of series."))
    end

    # Generates ad-hoc functions to compute the IAU Series including only 
    # the non-null coefficients. 

    # The elements of cpoly must always be >= ctrig
    fcn_body = Expr(:block)

    # Pre-computes all the Fundamental Arguments multiplications
    fad = Dict{Symbol,Vector{Int}}()
    for (m, fm) in enumerate(_FA_MAPPING)
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
    for j in 1:nseries
        var = Symbol("x$(j)x")

        # Automatically includes transformation to radians
        pcall[j] = Expr(
            :(=),
            var,
            Expr(
                :call,
                :(/),
                Expr(:call, :(*), Expr(:macrocall, Symbol("@evalpoly"), :, :t), π * 1e-6),
                648000,
            ),
        )
    end

    # Stores the indexes of all the contributions that have already been computed 
    dct = Dict{Tuple{Int,Int},Vector{Int}}()

    for j in 1:nseries
        # Assings all the polynomial contributions
        for i in eachindex(cpoly[j])
            dct[(j, i)] = Vector{Int}() # initialises dictionary keys 
            push!(fcn_body.args, Expr(:(=), Symbol("x$(j)x$i"), cpoly[j][i]))
            push!(pcall[j].args[2].args[2].args[2].args, Symbol("x$(j)x$i"))
        end
    end

    for jj in 1:nseries
        for i in eachindex(ctrig[jj])
            # Parses the trigonometric series
            for s in ctrig[jj][i]

                # Computes the ARGUMENT expression
                argv = Vector{Symbol}()
                for (j, n) in enumerate(s.N)
                    if n != 0
                        push!(
                            argv,
                            Symbol("$(_FA_MAPPING[j])$(argmin(abs.(fad[_FA_MAPPING[j]] .- n)))"),
                        )
                    end
                end

                arge = length(argv) > 0 ? Expr(:call, :(+), argv...) : 0.0

                # Stores current block expression
                blkₑ = Vector{Expr}()

                # Re-parses the entire trigonometric coefficients vector to check 
                # for elements with same arg coeffs 
                cn, sn = false, false

                for jjj in 1:nseries
                    for ii in eachindex(ctrig[jjj])
                        sc, cc = 0.0, 0.0

                        for (k, sk) in enumerate(ctrig[jjj][ii])
                            if sk.N == s.N && !(k in dct[(jjj, ii)])
                                # Stores the current index
                                push!(dct[(jjj, ii)], k)

                                # Groups together all the sin\cos terms
                                sc += sk.sc
                                cc += sk.cc

                                # Stores whether sin\cos are actually computed
                                sn = sn ? sn : sk.sc != 0.0
                                cn = cn ? cn : sk.cc != 0.0
                            end
                        end

                        sumₑ = Expr(:call, :(+))
                        sc != 0.0 && push!(sumₑ.args, Expr(:call, :(*), sc, :sarg))
                        cc != 0.0 && push!(sumₑ.args, Expr(:call, :(*), cc, :carg))

                        if sc != 0.0 || cc != 0.0
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

    for j in 1:nseries
        push!(fcn_body.args, pcall[j])
    end

    # Adds return types
    if nseries > 1
        rtupl = Expr(:tuple)
        for j in 1:nseries
            push!(rtupl.args, Symbol("x$(j)x"))
        end

        push!(fcn_body.args, Expr(:return, rtupl))
    end

    # Assembles the function
    fcn = _assemble_function(iau_model, fname, fcn_body)
    return eval(fcn)
end


function build_nutation_series(
    fname::Symbol, iau_model::Symbol, lon::AbstractVector{N}, obl::AbstractVector{N}
) where {N <: AbstractVector{<:IAUSeries}}

    # Generates ad-hoc functions to compute the IAU Series including only 
    # the non-null coefficients. 

    # ψseries and ϵseries must have the same number of elements (should be 2)
    nblk = length(lon)
    if nblk != length(obl)
        throw(ArgumentError("The longitude and obliquity series must have the same order."))
    end

    fbody = Expr(:block)

    # Pre-computes all the Fundamental Arguments multiplications
    fad = Dict{Symbol,Vector{Int}}()
    for (m, fm) in enumerate(_FA_MAPPING)

        fad[fm] = unique(vcat([map(x->x.N[m], s[i]) for i = 1:nblk for s in (lon, obl)]...))
        filter!(x-> x != 0, fad[fm])

        # Register the expression for that facotr
        for (k, mk) in enumerate(fad[fm])
            fae = Expr(:(.), :fa, QuoteNode(fm))
            push!(fbody.args, Expr(:(=), Symbol("$fm$k"), Expr(:call, :(*), mk, fae)))
        end
    end

    # Final @evalpoly call expression (depends on nblk)
    pcall = Vector{Expr}(undef, 2)
    for (j, var) in enumerate((:dψ, :dϵ))
        
        # Automatically includes transformation to radians
        pcall[j] = Expr(
            :(=),
            var,
            Expr(
                :call,
                :(/),
                Expr(:call, :(*), Expr(:macrocall, Symbol("@evalpoly"), :, :t), π * 1e-6),
                648000,
            ),
        )
    end

    # Stores the indexes of all the contributions that have already been computed 
    # First Tuple index refers to ψ or ϵ, the second to the corresponding block
    dct = Dict{Tuple{Int,Int},Vector{Int}}()

    # Initialises the variables 
    for i in 1:nblk
        # Adds ψ and ϵ contributions 
        for (j, var) in enumerate((Symbol("dψ$(i-1)"), Symbol("dϵ$(i-1)"))) 
            dct[(j, i)] = Vector{Int}() 
            push!(fbody.args, Expr(:(=), var, 0))
            push!(pcall[j].args[2].args[2].args[2].args, var)
        end
    end

    for i in 1:nblk
        for series in (lon, obl)

            # Parses ψ or ϵ block 
            for s in series[i]
                
                # for (j, n) in enumerate(s.N)
                #     if n != 0
                #         push!(
                #             argv,
                #             Symbol("$(_FA_MAPPING[j])$(argmin(abs.(fad[_FA_MAPPING[j]] .- n)))"),
                #         )
                #     end
                # end

                # Computes the ARGUMENT expression
                argv = Vector{Symbol}()
                for (fa, n) in zip(_FA_MAPPING, s.N)
                    if n != 0 
                        push!(argv, Symbol("$(fa)$(argmin(abs.(fad[fa] .- n)))"))
                    end
                end

                arge = length(argv) > 0 ? Expr(:call, :(+), argv...) : 0

                # Stores current block expression
                blkₑ = Vector{Expr}()

                # Re-parses the entire ψ and ϵ vectors to check 
                # for elements with same arg coeffs 
                cn, sn = false, false

                for ii in 1:nblk

                    # Parse ψ and ϵ
                    for (v, vct) in enumerate([lon, obl])
                        sc, cc = 0.0, 0.0

                        # Parses the inner series 
                        for (k, sk) in enumerate(vct[ii])
                            if sk.N == s.N && !(k in dct[(v, ii)])

                                # Stores the current index
                                push!(dct[(v, ii)], k)

                                # Groups together all the sin\cos terms
                                sc += sk.sc
                                cc += sk.cc

                                # Stores whether sin\cos are actually computed
                                sn = sn ? sn : sk.sc != 0.0
                                cn = cn ? cn : sk.cc != 0.0
                            end
                        end

                        sumₑ = Expr(:call, :(+))
                        sc != 0.0 && push!(sumₑ.args, Expr(:call, :(*), sc, :sarg))
                        cc != 0.0 && push!(sumₑ.args, Expr(:call, :(*), cc, :carg))

                        if sc != 0.0 || cc != 0.0
                            var = v == 1 ? "ψ" : "ϵ"
                            push!(blkₑ, Expr(:(+=), Symbol("d" * var * "$(ii-1)"), sumₑ))
                        end
                    end
                end

                # Adds all the expressions 
                if cn || sn

                    # Reassigns the arg variable
                    push!(fbody.args, Expr(:(=), :arg, arge))

                    # Computes and stores sin(ARG) and cos(ARG)
                    cn && push!(fbody.args, Expr(:(=), :carg, Expr(:call, :cos, :arg)))
                    sn && push!(fbody.args, Expr(:(=), :sarg, Expr(:call, :sin, :arg)))

                    push!(fbody.args, blkₑ...)
                end
            end
        end
    end

    for i in 1:2
        push!(fbody.args, pcall[i])
    end

    # Adds return types 
    push!(fbody.args, Expr(:return, Expr(:tuple, :dψ, :dϵ)))

    # Assembles the function
    fcn = _assemble_function(iau_model, fname, fbody)
    return eval(fcn)

end