
# This function is used to automatically generated optimised functions that 
# compute the series associated to the IAU constants.

function build_series(
    fname::Symbol, iau_model::Symbol, trig::AbstractVector, poly=nothing
)

    # Used to map the IAUSeries coefficients to the fundamental argument name 
    FA_MAPS = SVector(
        :Mₐ, :Sₐ, :uₘ, :Dₛ, :Ωₘ, :λ_Me, :λ_Ve, :λ_Ea, :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ
    )

    # Copy the expressions since it will be manipulated 
    ctrig = deepcopy(trig)

    # Retrieve number of series, i.e., the number of returned outputs
    nout = length(ctrig)
    if isnothing(poly)
        poly = [zeros(length(ctrig[i])) for i = 1:nout]
    elseif length(poly) != nout
        throw(ArgumentError("`trig` and `poly` must have the same length."))
    end

    # The elements of cpoly must always be >= ctrig
    fbody = Expr(:block)

    # Retrieve unique sets of coefficients multiplying the Fundamental Arguments 
    fa_list::Vector{SVector{14, Int}} = map(x->x.N, unique(x->x.N, vcat(vcat(ctrig...)...)))

    # Pre-computes all the Fundamental Arguments multiplications
    fad = Dict{Symbol, Vector{Int}}()
    for (m, fm) in enumerate(FA_MAPS)

        # Retrieve the unique coefficients multiplying the m-th fundamental argument.
        fad[fm] = filter(x-> x != 0, unique(map(x->x[m], fa_list)))

        # Skip the remaining if the argument is never used.
        len = length(fad[fm])
        len == 0 && continue

        # Register the expressions for that factor
        name = Expr(:(.), :fa, QuoteNode(fm))
        mult::Vector{Expr} = [
            Expr(:(=), Symbol("$fm$k"), Expr(:call, :*, fad[fm][k], name)) for k = 1:len
        ]

        append!(fbody.args, mult)

    end

    # Final @evalpoly call expression (depends on number of outputs)
    # Automatically includes transformation from arcseconds to radians
    epoly = Expr(:call, :*, Expr(:macrocall, Symbol("@evalpoly"), :, :t), 1e-6*π/648000)
    pcall = [Expr(:(=), Symbol("x$j"), deepcopy(epoly)) for j = 1:nout]

    # Assings all the polynomial contributions
    for j in 1:nout
        for i = 1:max(length(ctrig[j]), length(poly[j]))
            pval = i > length(poly[j]) ? 0.0 : poly[j][i]
            push!(fbody.args, Expr(:(=), Symbol("x$(j)x$i"), pval))
            push!(pcall[j].args[2].args[2].args, Symbol("x$(j)x$i"))
        end
    end

    # For each unique ARGUMENT expression 
    for fa_set in reverse(fa_list) 

        # Stores the expressions for this FA set
        setₑ = Vector{Expr}()

        cₙ, sₙ = false, false
        for (j, series) in enumerate(ctrig) 
            for (i, block) in enumerate(series)

                # There should always be single one! 
                idx = findfirst(x->x.N == fa_set, block)
                isnothing(idx) && continue

                # Remove this element to make future searches faster
                row = popat!(ctrig[j][i], idx)

                if row.sc != 0 || row.cc != 0 
                    sₙ = sₙ ? sₙ : row.sc != 0
                    cₙ = cₙ ? cₙ : row.cc != 0

                    sumₑ = Expr(:call, :(+))
                    row.sc != 0 && push!(sumₑ.args, Expr(:call, :*, row.sc, :sarg))
                    row.cc != 0 && push!(sumₑ.args, Expr(:call, :*, row.cc, :carg))
                    
                    push!(setₑ, Expr(:(+=), Symbol("x$(j)x$i"), sumₑ))
                end

            end
        end

        # Computes the sin\cos for the ARGUMENT expression
        if cₙ || sₙ
            # Retrieve the list of non-null FA  
            arg_list = Symbol[]
            for (j, fa) in enumerate(fa_set)
                if fa != 0
                    push!(
                        arg_list, 
                        Symbol("$(FA_MAPS[j])$(argmin(abs.(fad[FA_MAPS[j]] .- fa)))")
                    )
                end
            end

            # Reassigns the ARG variable 
            push!(fbody.args, Expr(:(=), :arg,  Expr(:call, :(+), arg_list...)))

            # Compute the sin(ARG) and cos(ARG)
            sₙ && push!(fbody.args, Expr(:(=), :sarg, Expr(:call, :sin, :arg)))
            cₙ && push!(fbody.args, Expr(:(=), :carg, Expr(:call, :cos, :arg)))

            push!(fbody.args, setₑ...)
        end
    end

    for j in 1:nout
        push!(fbody.args, pcall[j])
    end

    # Adds return types
    if nout > 1
        rtupl = Expr(:tuple)
        for j in 1:nout
            push!(rtupl.args, Symbol("x$j"))
        end

        push!(fbody.args, Expr(:return, rtupl))
    end

    # Assembles the function
    fcn = _assemble_function(iau_model, fname, fbody)
    return eval(fcn)
end

function _assemble_function(iau_model, fname, fbody)

    # Creates function call Expression
    fcall = Expr(:call, fname)

    push!(fcall.args, Expr(:(::), iau_model))   # Adds model
    push!(fcall.args, Expr(:(::), :t, :Number)) # Adds time 
    push!(fcall.args, Expr(:(::), :fa, :FundamentalArguments)) # Adds FA

    # Assembles Function head and body
    return Expr(:function, fcall, fbody)

end