# Reload natural function of NaturalSort package
function NaturalSort.natural(x::Real,y::Real)
    isless(x,y)
end

# Simplified implementation of 'smooth.splint' from R language.
function SmoothSpline(y::Vector{<: Real},
        df::Integer,
        x::Vector{<: Real} = nothing;
        w::Vector{<: Real} = nothing,
        tol::AbstractFloat = 1e-4,
        spar::AbstractFloat = NaN)

    if !isnothing(x)
        @assert length(x) == length(y)
    end
    if isnothing(w)
        w = repeat([1.],length(x))
    else
        @assert length(w) == length(x)
        @assert all(w .>= 0)
        @assert !all(w .== 0)
        w = (w .* sum(w .> 0)) ./ sum(w)
    end

    xx = round.((x .- mean(x)) ./ tol)
    iOx = issorted(x) ? repeat([true],length(x)) : sortperm(x)
    xxs = xx[iOx]
    nd = append!([true],xxs[1:(end - 1)] .< xxs[2:end])
    ux = x[iOx][nd]
    nx = length(ux)
    @assert nx > 3
    if nx == length(x)
        ox = issorted(x) ? true : 
            (
                xx = copy(iOx);
                xxx = copy(iOx);
                resize!(xx,maximum(xx));
                xx[xxx] .= 1:length(xxx);
                xx
            )
        tmp = hcat(w,w .* y,w .* y .^ 2)[iOx,:]
    else
        ox = indexin(xx,xxs[nd]) |> α -> Vector{Union{Nothing,Float64}}(α)
        ox[isnothing.(ox)] .= NaN
        tmp = [ (1:length(x))[j] for j in [ ox .== i for i in unique(ox) ] ] |> 
            α -> all(w .== 1) ? 
                [ vcat(length(β),sum(y[β]),sum(y[β] .^ 2)) for β in α ] : 
                [ vcat(sum(w[β]),sum(w[β] .* y[β]),sum(w[β] .* y[β] .^ 2)) 
                 for β in α ]
        tmp = hcat(tmp...)'
        tmp_idx = sortperm(tmp[:,2])
        tmp = Float64.(tmp[tmp_idx,:])
    end
    wbar = tmp[:,1]
    ybar = tmp[:,2] ./ [ i > 0 ? i : 1 for i in wbar ]
    yssw = sum(tmp[:,3] .- wbar .* ybar .^ 2)
    # clean variables
    iOx,xx,xxs,nd,tmp = nothing,nothing,nothing,nothing,nothing

    # cv/CV is false, all_knots is false (R function parameters)
    r_ux = ux[end] - ux[1]
    xbar = (ux .- ux[1]) ./ r_ux
    nknots = NknotsSmspl(nx)
    knot = vcat(repeat([xbar[1]],3),
                xbar[trunc.(Int,[ i for i in range(1;length=nknots,stop=nx) ])],
                repeat([xbar[end]],3))
    nk = nknots + 2

    if isnan(spar)
        spar = 0.
        ispar = 0
    else
        ispar = 1
    end
    icrit = 1
    dofoff = 0
    if df > 1 && df <= nx
        icrit = 3
        dofoff = df
    end

    iparms = [icrit,ispar,500,false]
    coef = zeros(nk)
    sz = zeros(nx)
    lev = zeros(nx)
    scrtch = zeros(18 * nk + 1)
    crit = Ref{Real}(0)
    parms = [-1.5,1.5,1e-4,2e-8,-1]
    ier = Ref{Real}(0)

    RBart!(Ref(1),Ref(dofoff),xbar,ybar,wbar,Ref(yssw),Ref(nx),knot,Ref(nk),
           coef,sz,lev,crit,iparms,Ref(spar),parms,scrtch,Ref(4),Ref(1),ier)

    return sz
end

function NknotsSmspl(n)::Integer
    if n < 50
        return n
    else
        a1,a2,a3,a4 = log2.([50,100,140,200])
        if n < 200
            return 2 ^ (a1 + (a2 - a1) * (n - 50) / 150) |> trunc
        elseif n < 800
            return 2 ^ (a2 + (a3 - a2) * (n - 200) / 600) |> trunc
        elseif n < 3200
            return 2 ^ (a3 + (a4 - a3) * (n - 800) / 2400) |> trunc
        else
            return 200 + (n - 3200) ^ 0.2 |> trunc
        end
    end
end

function FSign(x::Real,y::Real)
    if isnan(x) || isnan(y)
        return x + y
    end
    return y >= 0 ? abs(x) : -abs(x)
end

function CreateSBart()
    # Static local variables
    c_gold = 1 - 1 / MathConstants.golden
    ratio = 1.

    function InnerSBart!(penalt::Base.RefValue{<: Real},
            dofoff::Base.RefValue{<: Real},xs::Vector{<: Real},
            ys::Vector{<: Real},ws::Vector{<: Real},ssw::Base.RefValue{<: Real},
            n::Base.RefValue{<: Integer},knot::Vector{<: Real},
            nk::Base.RefValue{<: Integer},coef::Vector{<: Real},
            sz::Vector{<: Real},lev::Vector{<: Real},
            crit::Base.RefValue{<: Real},iparms::Vector{<: Integer},
            iparms_index::Vector{<: Integer},spar::Base.RefValue{<: Real},
            parms::Vector{<: Real},parms_index::Vector{<: Integer},
            isetup::Base.RefValue{<: Integer},scrtch::Vector{<: Real},
            scrtch_index::Vector{<: Integer},ld4::Base.RefValue{<: Integer},
            ldnk::Base.RefValue{<: Integer},ier::Base.RefValue{<: Real})

        icrit,ispar,iter = iparms_index
        lspar,uspar,tol,ϵ,Ratio = parms_index
        xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3 = scrtch_index[1:9]

        spar_is_λ = false
        d = fu = u = 0.

        ws[ws .> 0] .= sqrt.(ws[ws .> 0])
        if isetup[] < 0
            spar_is_λ = true
        elseif isetup[] != 1
            SGRAM!(scrtch,[sg0,sg1,sg2,sg3],knot,nk)
            STXWX!(xs,ys,ws,n,knot,nk,scrtch,[xwy,hs0,hs1,hs2,hs3])
            spar_is_λ = isetup[] == 2
            if !spar_is_λ
                t1 = sum(scrtch[(hs0 + 3 - 1):(hs0 + (nk[] - 3) - 1)])
                t2 = sum(scrtch[(sg0 + 3 - 1):(sg0 + (nk[] - 3) - 1)])
                ratio = t1 / t2
            end
            isetup[] = 1
        end

        if iparms[ispar] == 1
            parms[lspar] = spar_is_λ ? spar[] : (ratio * 16 ^ (spar[] * 6 - 2))
            SSLVRG!(penalt,dofoff,xs,ys,ws,ssw,n,knot,nk,coef,sz,lev,crit,
                    iparms[icrit],parms[lspar],scrtch,scrtch_index,ld4,ldnk,ier)
            parms[Ratio] = ratio
            return
        end

        ax = parms[lspar]
        bx = parms[uspar]
        maxit = iparms[iter]
        iparms[iter] = 0
        a = ax
        b = bx
        v = a + c_gold * (b - a)
        w = v
        x = v
        e = 0.
        parms[lspar] = spar_is_λ ? x : (ratio * 16 ^ (x * 6 - 2))
        SSLVRG!(penalt,dofoff,xs,ys,ws,ssw,n,knot,nk,coef,sz,lev,crit,
                iparms[icrit],parms[lspar],scrtch,scrtch_index,ld4,ldnk,ier)
        fx = crit[]
        fv = fx
        fw = fx

        # Ultimate loop
        while ier[] == 0
            xm = (a + b) / 2
            tol1 = parms[ϵ] * abs(x) + parms[tol] / 3
            tol2 = tol1 * 2
            iparms[iter] += 1
            if abs(x - xm) <= (tol2 - (b - a) * 0.5) || iparms[iter] > maxit
                # L_End
                parms[Ratio] = ratio
                spar[] = x
                crit[] = fx
                break
            end

            if abs(e) <= tol1 || isinf(fx) || isinf(fv) || isinf(fw)
                # L_GoldenSect
                if x >= xm
                    e = a - x
                else
                    e = b - x
                end
                d = c_gold * e
                # L50
                u = x + (abs(d) >= tol1 ? d : FSign(tol1,d))
                parms[lspar] = spar_is_λ ? u : (ratio * 16 ^ (u * 6 - 2))
                SSLVRG!(penalt,dofoff,xs,ys,ws,ssw,n,knot,nk,coef,sz,lev,crit,
                        iparms[icrit],parms[lspar],scrtch,scrtch_index,ld4,ldnk,
                        ier)
                fu = crit[]
                if fu <= fx
                    if u >= x
                        a = x
                    else
                        b = x
                    end
                    v,fv = w,fw
                    w,fw = x,fx
                    x,fx = u,fu
                else
                    if u < x
                        a = u
                    else
                        b = u
                    end
                    if fu <= fw || w == x
                        v,fv = w,fw
                        w,fw = u,fu
                    elseif fu <= fv || v == x || v == w
                        v,fv = u,fu
                    end
                end
                continue
            end

            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2
            if q > 0
                p = -p
            end
            q = abs(q)
            r = e
            e = d
            if abs(p) >= abs(0.5 * q * r) || q == 0
                # L_GoldenSect
                if x >= xm
                    e = a - x
                else
                    e = b - x
                end
                d = c_gold * e
                # L50
                u = x + (abs(d) >= tol1 ? d : FSign(tol1,d))
                parms[lspar] = spar_is_λ ? u : (ratio * 16 ^ (u * 6 - 2))
                SSLVRG!(penalt,dofoff,xs,ys,ws,ssw,n,knot,nk,coef,sz,lev,crit,
                        iparms[icrit],parms[lspar],scrtch,scrtch_index,ld4,ldnk,
                        ier)
                fu = crit[]
                if fu <= fx
                    if u >= x
                        a = x
                    else
                        b = x
                    end
                    v,fv = w,fw
                    w,fw = x,fx
                    x,fx = u,fu
                else
                    if u < x
                        a = u
                    else
                        b = u
                    end
                    if fu <= fw || w == x
                        v,fv = w,fw
                        w,fw = u,fu
                    elseif fu <= fv || v == x || v == w
                        v,fv = u,fu
                    end
                end
                continue
            end
            
            if p <= q * (a - x) || p >= q * (b - x)
                # L_GoldenSect
                if x >= xm
                    e = a - x
                else
                    e = b - x
                end
                d = c_gold * e
                # L50
                u = x + (abs(d) >= tol1 ? d : FSign(tol1,d))
                parms[lspar] = spar_is_λ ? u : (ratio * 16 ^ (u * 6 - 2))
                SSLVRG!(penalt,dofoff,xs,ys,ws,ssw,n,knot,nk,coef,sz,lev,crit,
                        iparms[icrit],parms[lspar],scrtch,scrtch_index,ld4,ldnk,
                        ier)
                fu = crit[]
                if fu <= fx
                    if u >= x
                        a = x
                    else
                        b = x
                    end
                    v,fv = w,fw
                    w,fw = x,fx
                    x,fx = u,fu
                else
                    if u < x
                        a = u
                    else
                        b = u
                    end
                    if fu <= fw || w == x
                        v,fv = w,fw
                        w,fw = u,fu
                    elseif fu <= fv || v == x || v == w
                        v,fv = u,fu
                    end
                end
                continue
            end

            d = p / q
            u = x + d
            if (u - a) < tol2 || (b - u) < tol2
                d = FSign(tol1,xm - x)
            end
            # L50
            u = x + (abs(d) >= tol1 ? d : FSign(tol1,d))
            parms[lspar] = spar_is_λ ? u : (ratio * 16 ^ (u * 6 - 2))
            SSLVRG!(penalt,dofoff,xs,ys,ws,ssw,n,knot,nk,coef,sz,lev,crit,
                    iparms[icrit],parms[lspar],scrtch,scrtch_index,ld4,ldnk,
                    ier)
            fu = crit[]
            if fu <= fx
                if u >= x
                    a = x
                else
                    b = x
                end
                v,fv = w,fw
                w,fw = x,fx
                x,fx = u,fu
            else
                if u < x
                    a = u
                else
                    b = u
                end
                if fu <= fw || w == x
                    v,fv = w,fw
                    w,fw = u,fu
                elseif fu <= fv || v == x || v == w
                    v,fv = u,fu
                end
            end
        end
    end # End of function
end

SBart! = CreateSBart()

function RBart!(penalt::Base.RefValue{<: Real},dofoff::Base.RefValue{<: Real},
        xs::Vector{<: Real},ys::Vector{<: Real},ws::Vector{<: Real},
        ssw::Base.RefValue{<: Real},n::Base.RefValue{<: Integer},
        knot::Vector{<: Real},nk::Base.RefValue{<: Integer},
        coef::Vector{<: Real},sz::Vector{<: Real},lev::Vector{<: Real},
        crit::Base.RefValue{<: Real},iparms::Vector{<: Integer},
        spar::Base.RefValue{<: Real},parms::Vector{<: Real},
        scrtch::Vector{<: Real},ld4::Base.RefValue{<: Integer},
        ldnk::Base.RefValue{<: Integer},ier::Base.RefValue{<: Real})

    @assert length(iparms) == 4
    @assert length(xs) == n[]
    @assert length(ys) == n[]
    @assert length(ws) == n[]
    @assert length(knot) == (nk[] + 4)
    @assert length(coef) == nk[]
    @assert length(sz) == n[]
    @assert length(lev) == n[]
    @assert length(parms) == 5

    if iparms[4] == 1
        isetup = Ref(2)
    else
        isetup = Ref(0)
    end

    SBart!(penalt,dofoff,xs,ys,ws,ssw,n,knot,nk,coef,sz,lev,crit,iparms,[1,2,3],
           spar,parms,[1,2,3,4,5],isetup,scrtch,
           [1,nk[] + 1,2 * nk[] + 1,3 * nk[] + 1,4 * nk[] + 1,5 * nk[] + 1,
            6 * nk[] + 1,7 * nk[] + 1,8 * nk[] + 1,9 * nk[] + 1,
            9 * nk[] + ld4[] * nk[] + 1,9 * nk[] + 2 * ld4[] * nk[] + 1],
           ld4,ldnk,ier)
end

function SGRAM!(scrtch::Vector{<: Real},
        scrtch_index::Vector{<: Integer},
        tb::Vector{<: Real},
        nb::Base.RefValue{<: Integer})

    @assert length(tb) == (nb[] + 4)
    sg0,sg1,sg2,sg3 = scrtch_index

    mflag = Ref(0)
    work = Vector{Float64}(undef,16)
    vnikx = Matrix{Float64}(undef,4,3)

    scrtch[sg0:(sg0 + nb[] - 1)] .= scrtch[sg1:(sg1 + nb[] - 1)] .= 0
    scrtch[sg2:(sg2 + nb[] - 1)] .= scrtch[sg3:(sg3 + nb[] - 1)] .= 0
    ileft = 1

    for i in 1:nb[]
        ileft = FindInterval!(tb,nb[] + 1,tb[i],false,false,ileft,mflag)
        BSPLVD!(tb,4,tb[i],ileft,work,vnikx,3)
        yw1 = [ vnikx[i,3] for i in 1:4 ]
        BSPLVD!(tb,4,tb[i + 1],ileft,work,vnikx,3)
        yw2 = [ vnikx[i,3] - yw1[i] for i in 1:4 ]
        wpt = tb[i + 1] - tb[i]

        if ileft >= 4
            for ii in 1:4
                scrtch[sg0 + (ileft - 4 + ii) - 1] += wpt * 
                    (yw1[ii] ^ 2 + 
                     (yw2[ii] * yw1[ii] + yw2[ii] * yw1[ii]) * 0.5 + 
                     yw2[ii] ^ 2 * 0.333)
                jj = ii + 1
                if jj <= 4
                    scrtch[sg1 + (ileft + ii - 4) - 1] += wpt * 
                        (yw1[ii] * yw1[jj] + 
                         (yw2[ii] * yw1[jj] + yw2[jj] * yw1[ii]) * 0.5 + 
                         yw2[ii] * yw2[jj] * 0.333)
                end
                jj += 1
                if jj <= 4
                    scrtch[sg2 + (ileft + ii - 4) - 1] += wpt * 
                        (yw1[ii] * yw1[jj] + 
                         (yw2[ii] * yw1[jj] + yw2[jj] * yw1[ii]) * 0.5 + 
                         yw2[ii] * yw2[jj] * 0.333)
                end
                jj += 1
                if jj <= 4
                    scrtch[sg3 + (ileft + ii - 4) - 1] += wpt * 
                        (yw1[ii] * yw1[jj] + 
                         (yw2[ii] * yw1[jj] + yw2[jj] * yw1[ii]) * 0.5 + 
                         yw2[ii] * yw2[jj] * 0.333)
                end
            end
        elseif ileft == 3
            for ii in 1:3
                scrtch[sg0 + (ileft - 3 + ii) - 1] += wpt * 
                    (yw1[ii] ^ 2 + 
                     (yw2[ii] * yw1[ii] + yw2[ii] * yw1[ii]) * 0.5 + 
                     yw2[ii] ^ 2 * 0.333)
                jj = ii + 1
                if jj <= 3
                    scrtch[sg1 + (ileft + ii - 3) - 1] += wpt * 
                        (yw1[ii] * yw1[jj] + 
                         (yw2[ii] * yw1[jj] + yw2[jj] * yw1[ii]) * 0.5 + 
                         yw2[ii] * yw2[jj] * 0.333)
                end
                jj += 1
                if jj <= 3
                    scrtch[sg2 + (ileft + ii - 3) - 1] += wpt * 
                    (yw1[ii] * yw1[jj] + 
                     (yw2[ii] * yw1[jj] + yw2[jj] * yw1[ii]) * 0.5 + 
                     yw2[ii] * yw2[jj] * 0.333)
                end
            end
        elseif ileft == 2
            for ii in 1:2
                scrtch[sg0 + (ileft - 2 + ii) - 1] += wpt * 
                    (yw1[ii] ^ 2 + 
                     (yw2[ii] * yw1[ii] + yw2[ii] * yw1[ii]) * 0.5 + 
                     yw2[ii] ^ 2 * 0.333)
                jj = ii + 1
                if jj <= 2
                    scrtch[sg1 + (ileft + ii - 2) - 1] += wpt * 
                        (yw1[ii] * yw1[jj] + 
                         (yw2[ii] * yw1[jj] + yw2[jj] * yw1[ii]) * 0.5 + 
                         yw2[ii] * yw2[jj] * 0.333)
                end
            end
        elseif ileft == 1
            ii = jj = 1
            scrtch[sg0 + (ileft - 1 + ii) - 1] += wpt * 
                (yw1[ii] ^ 2 + 
                 (yw2[ii] * yw1[ii] + yw2[ii] * yw1[ii]) * 0.5 + 
                 yw2[ii] ^ 2 * 0.333)
        end
    end
end

function STXWX!(x::Vector{<: Real},z::Vector{<: Real},w::Vector{<: Real},
        k::Base.RefValue{<: Integer},xknot::Vector{<: Real},
        n::Base.RefValue{<: Integer},scrtch::Vector{<: Real},
        scrtch_index::Vector{<: Integer})

    @assert length(x) == k[]
    @assert length(z) == k[]
    @assert length(w) == k[]
    @assert length(xknot) == (n[] + 4)

    y,hs0,hs1,hs2,hs3 = scrtch_index

    scrtch[y:(y + n[] - 1)] .= scrtch[hs0:(hs0 + n[] - 1)] .= 0
    scrtch[hs1:(hs1 + n[] - 1)] .= scrtch[hs2:(hs2 + n[] - 1)] .= 0
    scrtch[hs3:(hs3 + n[] - 1)] .= 0
    ileft = 1
    ϵ = 0.1e-9

    mflag = Ref(0)
    work = Vector{Float64}(undef,16)
    vnikx = Matrix{Float64}(undef,4,1)

    for i in 1:k[]
        ileft = FindInterval!(xknot,n[] + 1,x[i],false,false,ileft,mflag)
        if mflag[] == 1
            if x[i] <= (xknot[ileft] + ϵ)
                ileft -= 1
            else
                return
            end
        end
        BSPLVD!(xknot,4,x[i],ileft,work,vnikx,1)
        j = ileft - 3
        scrtch[y + j - 1] += w[i] ^ 2 * z[i] * vnikx[1,1]
        scrtch[hs0 + j - 1] += w[i] ^ 2 * vnikx[1,1] ^ 2
        scrtch[hs1 + j - 1] += w[i] ^ 2 * vnikx[1,1] * vnikx[2,1]
        scrtch[hs2 + j - 1] += w[i] ^ 2 * vnikx[1,1] * vnikx[3,1]
        scrtch[hs3 + j - 1] += w[i] ^ 2 * vnikx[1,1] * vnikx[4,1]
        j = ileft - 2
        scrtch[y + j - 1] += w[i] ^ 2 * z[i] * vnikx[2,1]
        scrtch[hs0 + j - 1] += w[i] ^ 2 * vnikx[2,1] ^ 2
        scrtch[hs1 + j - 1] += w[i] ^ 2 * vnikx[2,1] * vnikx[3,1]
        scrtch[hs2 + j - 1] += w[i] ^ 2 * vnikx[2,1] * vnikx[4,1]
        j = ileft - 1
        scrtch[y + j - 1] += w[i] ^ 2 * z[i] * vnikx[3,1]
        scrtch[hs0 + j - 1] += w[i] ^ 2 * vnikx[3,1] ^ 2
        scrtch[hs1 + j - 1] += w[i] ^ 2 * vnikx[3,1] * vnikx[4,1]
        j = ileft
        scrtch[y + j - 1] += w[i] ^ 2 * z[i] * vnikx[4,1]
        scrtch[hs0 + j - 1] += w[i] ^ 2 * vnikx[4,1] ^ 2
    end
end

function SSLVRG!(penalt::Base.RefValue{<: Real},dofoff::Base.RefValue{<: Real},
        x::Vector{<: Real},y::Vector{<: Real},w::Vector{<: Real},
        ssw::Base.RefValue{<: Real},n::Base.RefValue{<: Integer},
        knot::Vector{<: Real},nk::Base.RefValue{<: Integer},
        coef::Vector{<: Real},sz::Vector{<: Real},lev::Vector{<: Real},
        crit::Base.RefValue{<: Real},icrit::Integer,λ::Real,
        scrtch::Vector{<: Real},scrtch_index::Vector{<: Integer},
        ld4::Base.RefValue{<: Integer},ldnk::Base.RefValue{<: Integer},
        info::Base.RefValue{<: Real})

    @assert length(x) == n[]
    @assert length(y) == n[]
    @assert length(w) == n[]
    @assert length(knot) == nk[] + 4
    @assert length(coef) == nk[]
    @assert length(sz) == n[]
    @assert length(lev) == n[]

    xwy,hs0,hs1,hs2,hs3,sg0,sg1,sg2,sg3 = scrtch_index[1:9]

    abd,p1ip = scrtch_index[10:11]
    abd_p1ip_dict = [ (i,j) for j in 1:nk[] for i in 1:ld4[] ] |> 
        x -> Dict(x .=> eachindex(x))
    p2ip = scrtch_index[12]

    work = Vector{Float64}(undef,16)
    vnikx = Matrix{Float64}(undef,4,1)
    
    mflag = Ref(0)
    ileft = 1
    ϵ = 0.1 ^ 11
    
    coef[1:nk[]] .= scrtch[xwy:(xwy + nk[] - 1)]

    scrtch[abd .+ [ abd_p1ip_dict[(4,i)] for i in 1:nk[] ] .- 1] .= 
        scrtch[hs0 .+ collect(1:nk[]) .- 1] .+ 
        λ .* scrtch[sg0 .+ collect(1:nk[]) .- 1]
    scrtch[abd .+ [ abd_p1ip_dict[(3,i)] for i in 2:nk[] ] .- 1] .= 
        scrtch[hs1 .+ collect(1:(nk[] - 1)) .- 1] .+ 
        λ .* scrtch[sg1 .+ collect(1:(nk[] - 1)) .- 1]
    scrtch[abd .+ [ abd_p1ip_dict[(2,i)] for i in 3:nk[] ] .- 1] .= 
        scrtch[hs2 .+ collect(1:(nk[] - 2)) .- 1] .+ 
        λ .* scrtch[sg2 .+ collect(1:(nk[] - 2)) .- 1]
    scrtch[abd .+ [ abd_p1ip_dict[(1,i)] for i in 4:nk[] ] .- 1] .= 
        scrtch[hs3 .+ collect(1:(nk[] - 3)) .- 1] .+ 
        λ .* scrtch[sg3 .+ collect(1:(nk[] - 3)) .- 1]
    DPBFA!(scrtch,abd,ld4,nk,3,info)

    if info[] != 0
        return
    end

    DPBSL!(scrtch,abd,ld4,nk,3,coef)
    sz[1:n[]] .= [ BValue(knot,coef,nk[],4,x[i],0) for i in 1:n[] ]
    if icrit >= 1
        SINERP!(scrtch,abd,ld4,nk,p1ip,p2ip,ldnk,0)
        for i in 1:n[]
            xv = x[i]
            ileft = FindInterval!(knot,nk[] + 1,xv,false,false,ileft,mflag)
            if mflag[] == -1
                ileft = 4
                xv = knot[4] + ϵ
            elseif mflag[] == 1
                ileft = nk[]
                xv = knot[nk[] + 1] - ϵ
            end
            j = ileft - 3
            BSPLVD!(knot,4,xv,ileft,work,vnikx,1)
            b0,b1,b2,b3 = vnikx[1:4]
            lev[i] = sum([1,2,2,2,1,2,2,1,2,1] .* 
                         scrtch[p1ip .+ [ abd_p1ip_dict[k] 
                                         for k in 
                                         diag([ (i,j) 
                                               for i in 
                                               vcat(4:-1:1,4:-1:2,4:-1:3,4),
                                               j in 
                                               vcat([ repeat([j + i - 1],
                                                             4 - i + 1) 
                                                     for i in 1:4 ]...) ]) ] .- 
                                1] .* 
                         vcat([repeat([b0],4),repeat([b1],3),repeat([b2],2),
                               b3]...) .* 
                         [b0,b1,b2,b3,b1,b2,b3,b2,b3,b3]) * w[i] ^ 2
        end

        df = 0
        if icrit == 1
            rss = ssw[]
            sumw = 0
            rss += sum(((y[1:n[]] .- sz[1:n[]]) .* w[1:n[]]) .^ 2)
            df += sum(lev[1:n[]])
            sumw += sum(w[1:n[]] .^ 2)
            crit[] = (rss / sumw) / (1 - (dofoff[] + penalt[] * df) / sumw) ^ 2
        elseif icrit == 2
            crit[] = 0
            crit[] += sum(((y[1:n[]] .- sz[1:n[]]) .* w[1:n[]] ./ 
                           (1 .- lev[1:n[]])) .^ 2)
            crit[] /= n[]
        else
            df += sum(lev[1:n[]])
            if icrit == 3
                crit[] = 3 + (dofoff[] - df) ^ 2
            else
                crit[] = df - dofoff[]
            end
        end
    end
end

function FindInterval!(xt::Vector{<: Real},
        n::Integer,
        x::Real,
        rightmost_closed::Bool,
        all_inside::Bool,
        ilo::Integer,
        mflag::Base.RefValue{<: Integer};
        left_open::Bool=false)

    if n == 0
        mflag[] = 0
        return 0
    end

    if ilo <= 0
        if x < xt[1] || (left_open && x <= xt[1])
            # left_boundary
            mflag[] = -1
            return (all_inside || (rightmost_closed && x == xt[1])) ? 1 : 0
        end
        ilo = 1
    end

    ihi = ilo + 1
    if ihi >= n
        if x > xt[n] || (!left_open && x >= xt[n])
            # right_boundary
            mflag[] = 1
            return (all_inside || (rightmost_closed && x  == xt[n])) ? 
                (n - 1) : 
                n
        end
        if n <= 1
            mflag[] = -1
            return (all_inside || (rightmost_closed && x == xt[1])) ? 1 : 0
        end
        ilo = n - 1
        ihi = n
    end

    goto_label = nothing
    if x < xt[ihi] || (left_open && x <= xt[ihi])
        if x > xt[ilo] || (!left_open && x >= xt[ilo])
            mflag[] = 0
            return ilo
        end
        # decrease ilo to capture x
        if !left_open
            istep = 1
            while true
                ihi = ilo
                ilo -= istep
                if ilo <= 1
                    break
                end
                if x >= xt[ilo]
                    goto_label = 50
                    break
                end
                istep *= 2
            end
        else
            istep = 1
            while true
                ihi = ilo
                ilo -= istep
                if ilo <= 1
                    break
                end
                if x > xt[ilo]
                    goto_label = 51
                    break
                end
                istep *= 2
            end
        end
        ilo = 1
        if x < xt[1] || (left_open && x <= xt[1])
            mflag[] = -1
            return (all_inside || (rightmost_closed && x == xt[1])) ? 1 : 0
        end
    else
        # increase ihi to captrue x
        if !left_open
            istep = 1
            while true
                ilo = ihi
                ihi += istep
                if ihi >= n
                    break
                end
                if x < xt[ihi]
                    goto_label = 50
                    break
                end
                istep *= 2
            end
        else
            istep = 1
            while true
                ilo = ihi
                ihi += istep
                if ihi >= n
                    break
                end
                if x <= xt[ihi]
                    goto_label = 51
                    break
                end
                istep *= 2
            end
        end
        if x > xt[n] || (!left_open && x >= xt[n])
            mflag[] = 1
            return (all_inside || (rightmost_closed && x  == xt[n])) ? 
                (n - 1) : 
                n
        end
        ihi = n
    end
    
    if left_open
        goto_label = 51
    end
    
    if goto_label == 50
        while true
            middle = (ilo + ihi) ÷ 2
            if middle == ilo
                mflag[] = 0
                return ilo
            end
            if x >= xt[middle]
                ilo = middle
            else
                ihi = middle
            end
        end
    else # goto_label == 51
        while true
            middle = (ilo + ihi) ÷ 2
            if middle == ilo
                mflag[] = 0
                return ilo
            end
            if x > xt[middle]
                ilo = middle
            else
                ihi = middle
            end
        end
    end
end

function CreateBSPLVB()
    # Static local variables
    j = 1
    δ_l = Vector{Float64}(undef,20)
    δ_r = Vector{Float64}(undef,20)

    function InnerBSPLVB!(t::Vector{<: Real},jhigh::Integer,index::Integer,
            x::Real,left::Integer,biatx::Vector{<: Real})
        i = jp1 = 0
        saved = term = 0
        
        if index != 2
            j = 1
            biatx[1] = 1
            if j >= jhigh
                return
            end
        end
        
        while true
            jp1 = j + 1
            δ_r[j] = t[left + j] - x
            δ_l[j] = x - t[left + 1 - j]
            saved = 0
            for i in 1:j
                term = biatx[i] / (δ_r[i] + δ_l[jp1 - i])
                biatx[i] = saved + δ_r[i] * term
                saved = δ_l[jp1 - i] * term
            end
            biatx[jp1] = saved
            j = jp1
            if j >= jhigh
                return
            end
        end
    end # End of function
end

BSPLVB! = CreateBSPLVB()

function BSPLVD!(t::Vector{<: Real},k::Integer,x::Real,left::Integer,
        a::Vector{<: Real},dbiatx::Matrix{<: Real},nderiv::Integer)

    @assert size(dbiatx) == (k,nderiv)
    @assert length(a) == k * k

    # We can also set the 'a' as a matrix parameter of the function, 
    # and the 'reshape...' below is not necessary!
    # We don't know why the original fortran code was a vector, not a matrix!
    # 'reshape...' can keep the references to the elements!
    b = reshape(a,k,k)

    mhigh = max(min(nderiv,k),1)
    kp1 = k + 1
    # 'vec...' also can keep the references to the elements!
    BSPLVB!(t,kp1 - mhigh,1,x,left,vec(dbiatx))
    if mhigh == 1
        return
    end
    ideriv = mhigh
    for _ in 2:mhigh
        jp1mid = 1
        for j in ideriv:k
            dbiatx[j,ideriv] = dbiatx[jp1mid,1]
            jp1mid += 1
        end
        ideriv -= 1
        BSPLVB!(t,kp1 - ideriv,2,x,left,vec(dbiatx))
    end

    jlow = 1
    for i in 1:k
        for j in jlow:k
            b[j,i] = 0
        end
        jlow = i
        b[i,i] = 1
    end
    for m in 2:mhigh
        kp1mm = kp1 - m
        il = left
        i = k
        for _ in 1:kp1mm
            factor = kp1mm / (t[il + kp1mm] - t[il])
            for j in 1:i
                b[i,j] = (b[i,j] - b[i - 1,j]) * factor
            end
            il -= 1
            i -= 1
        end
        for i in 1:k
            s = 0.
            jlow = max(i,m)
            for j in jlow:k
                s += b[j,i] * dbiatx[j,m]
            end
            dbiatx[i,m] = s
        end
    end
end

function DPBFA!(scrtch::Vector{<: Real},abd::Integer,
        lda::Base.RefValue{<: Integer},n::Base.RefValue{<: Integer},m::Integer,
        info::Base.RefValue{<: Real})

    abd_dict = [ (i,j) for j in 1:n[] for i in 1:lda[] ] |> 
        x -> Dict(x .=> eachindex(x))

    for j in 1:n[]
        info[] = j
        s = 0.
        ik = m + 1
        jk = max(j - m,1)
        μ = max(m + 2 - j,1)
        if m >= μ
            for k in μ:m
                idx_1 = abd_dict[(ik,jk)]
                idx_2 = abd_dict[(μ,j)]
                # 'dot' is from LinearAlgebra
                t = scrtch[abd + abd_dict[(k,j)] - 1] - 
                    dot(scrtch[abd .+ collect(idx_1:(idx_1 + k - μ - 1)) .- 1],
                        scrtch[abd .+ collect(idx_2:(idx_2 + k - μ - 1)) .- 1]) 
                t /= scrtch[abd + abd_dict[(m + 1,jk)] - 1]
                scrtch[abd + abd_dict[(k,j)] - 1] = t
                s += t ^ 2
                ik -= 1
                jk += 1
            end
        end
        
        s = scrtch[abd + abd_dict[(m + 1,j)] - 1] - s
        if s <= 0
            return
        end
        scrtch[abd + abd_dict[(m + 1,j)] - 1] = sqrt(s)
    end
    info[] = 0
end

function DPBSL!(scrtch::Vector{<: Real},abd::Integer,
        lda::Base.RefValue{<: Integer},n::Base.RefValue{<: Integer},m::Integer,
        b::Vector{<: Real})

    abd_dict = [ (i,j) for j in 1:n[] for i in 1:lda[] ] |> 
        x -> Dict(x .=> eachindex(x))
    for k in 1:n[]
        lm = min(k - 1,m)
        la = m + 1 - lm
        lb = k - lm
        idx = abd_dict[(la,k)]
        t = dot(scrtch[abd .+ collect(idx:(idx + lm - 1)) .- 1],
                b[lb:(lb + lm - 1)])
        b[k] = (b[k] - t) / scrtch[abd + abd_dict[(m + 1,k)] - 1]
    end
    for kb in 1:n[]
        k = n[] + 1 - kb
        lm = min(k - 1,m)
        la = m + 1 - lm
        lb = k - lm
        b[k] /= scrtch[abd + abd_dict[(m + 1,k)] - 1]
        t = -b[k]
        idx = abd_dict[(la,k)]
        # 'axpy!' can NOT change the subset of array in place.
        # So we can use the macro '@view'
        axpy!(t,scrtch[abd .+ collect(idx:(idx + lm - 1)) .- 1],
              @view(b[lb:(lb + lm - 1)]))
    end
end

function SINERP!(scrtch::Vector{<: Real},abd::Integer,
        ld4::Base.RefValue{<: Integer},nk::Base.RefValue{<: Integer},
        p1ip::Integer,p2ip::Integer,ldnk::Base.RefValue{<: Integer},
        flag::Integer)
    
    abd_p1ip_dict = [ (i,j) for j in 1:nk[] for i in 1:ld4[] ] |> 
        x -> Dict(x .=> eachindex(x))
    p2ip_dict = [ (i,j) for j in 1:nk[] for i in ldnk[] ] |> 
        x -> Dict(x .=> eachindex(x))
    
    cv = zeros(4)
    wjm1,wjm2,wjm3 = zeros(1),zeros(2),zeros(3)
    
    for i = 1:nk[]
        j = nk[] - i + 1
        cv[1] = 1 / scrtch[abd + abd_p1ip_dict[(4,j)] - 1]
        if j <= (nk[] - 3)
            cv[2:4] .= scrtch[abd .+ 
                              diag([ abd_p1ip_dict[(tmp1,j + tmp2)] 
                                    for tmp1 in 1:3,tmp2 in 3:-1:1 ]) .- 
                              1] .* cv[1]
        elseif j == (nk[] - 2)
            cv[2] = 0
            cv[3:4] .= scrtch[abd .+ 
                              diag([ abd_p1ip_dict[(tmp1,j + tmp2)] 
                                    for tmp1 in 2:3,tmp2 in 2:-1:1 ]) .- 
                              1] .* cv[1]
        elseif j == (nk[] - 1)
            cv[2] = cv[3] = 0
            cv[4] = scrtch[abd + abd_p1ip_dict[(3,j + 1)] - 1] * cv[1]
        elseif j == nk[]
            cv[2] = cv[3] = cv[4] = 0
        end
        scrtch[p1ip + abd_p1ip_dict[(1,j)] - 1] = 0 - sum(cv[2:4] .* wjm3)
        scrtch[p1ip + abd_p1ip_dict[(2,j)] - 1] = 
            0 - sum(cv[2:4] .* [wjm3[2],wjm2[1],wjm2[2]])
        scrtch[p1ip + abd_p1ip_dict[(3,j)] - 1] = 
            0 - sum(cv[2:4] .* [wjm3[3],wjm2[2],wjm1[1]])
        scrtch[p1ip + abd_p1ip_dict[(4,j)] - 1] = 
            cv[1] ^ 2 + cv[2] ^ 2 * wjm3[1] + 2 * cv[2] * cv[3] * wjm3[2] + 
            2 * cv[2] * cv[4] * wjm3[3] + cv[3] ^ 2 * wjm2[1] + 
            2 * cv[3] * cv[4] * wjm2[2] + cv[4] ^ 2 * wjm1[1]
        wjm3[1] = wjm2[1]
        wjm3[2] = wjm2[2]
        wjm3[3] = scrtch[p1ip + abd_p1ip_dict[(2,j)] - 1]
        wjm2[1] = wjm1[1]
        wjm2[2] = scrtch[p1ip + abd_p1ip_dict[(3,j)] - 1]
        wjm1[1] = scrtch[p1ip + abd_p1ip_dict[(4,j)] - 1]
    end
    
    if flag != 0
        for i in 1:nk[]
            j = nk[] - i + 1
            for k in 1:4
                if (j + k - 1) <= nk[]
                    scrtch[p2ip + p2ip_dict[(j,j + k - 1)] - 1] = 
                        scrtch[p1ip + abd_p1ip_dict[(5 - k,j)] - 1]
                end
            end
        end
        for i in 1:nk[]
            j = nk[] - i + 1
            if (j - 4) >= 1
                for k in (j - 4):-1:1
                    cv[1] = 1 / scrtch[abd + abd_p1ip_dict[(4,k)] - 1]
                    cv[2:4] .= scrtch[abd .+ 
                                      diag([ abd_p1ip_dict[(tmp1,k + tmp2)] 
                                            for tmp1 in 1:3,tmp2 in 3:-1:1 ]) .- 
                                      1] .* cv[1]
                    idx = p2ip_dict[(k + 3,j)]
                    scrtch[p2ip + p2ip_dict[(k,j)] - 1] = 
                        0 - sum(cv[2:4] .* 
                                scrtch[p2ip .+ 
                                       [ i for i in idx:-1:(idx - 2) ] .- 1])
                end
            end
        end
    end
end

function BValue(t::Vector{<: Real},bcoef::Vector{<: Real},n::Integer,k::Integer,
        x::Real,jderiv::Integer)
    @assert length(bcoef) == n
    
    # Max degree of B-spline is 20
    aj = Vector{Float64}(undef,20)
    dm = Vector{Float64}(undef,20)
    dp = Vector{Float64}(undef,20)
    # Error flag
    mflag = Ref(0)

    # Check input derivative less than the degree of B-spline
    if jderiv >= k
        return 0
    end

    i = 1
    if x != t[n + 1] || t[n + 1] != t[n + k]
        i = FindInterval!(t,n + k,x,false,false,i,mflag)
        if mflag[] != 0
            @warn "BValue: 'mflag' should be 0!"
            return 0
        end
    else
        i = n
    end
    
    km1 = k - 1
    if km1 <= 0
        return bcoef[i]
    end
    
    jcmin = 1
    imk = i - k
    if imk >= 0
        dm[1:km1] .= x .- t[i .+ 1 .- collect(1:km1)]
    else
        jcmin = 1 - imk
        dm[1:i] .= x .- t[i .+ 1 .- collect(1:i)]
        aj[k .- collect(i:km1)] .= 0
        dm[i:km1] .= dm[i]
    end
    
    jcmax = k
    nmi = n - i
    if nmi >= 0
        dp[1:km1] .= t[i .+ collect(1:km1)] .- x
    else
        jcmax = k + nmi
        dp[1:jcmax] .= t[i .+ collect(1:jcmax)] .- x
        aj[collect(jcmax:km1) .+ 1] .= 0
        dp[(jcmax + 1):km1] .= dp[jcmax]
    end
    
    aj[jcmin:jcmax] .= bcoef[imk .+ collect(jcmin:jcmax)]
    
   if jderiv >= 1
        for j in 1:jderiv
            kmj = k - j
            aj[1:kmj] .= (aj[2:(kmj + 1)] .- aj[1:kmj]) ./ 
                (dm[kmj:-1:1] .+ dp[1:kmj]) .* 
                kmj
        end
    end
    
    if jderiv != km1
        for j in (jderiv + 1):km1
            kmj = k - j
            aj[1:kmj] .= (aj[2:(kmj + 1)] .* dm[kmj:-1:1] .+ 
                          aj[1:kmj] .* dp[1:kmj]) ./ 
                (dm[kmj:-1:1] .+ dp[1:kmj])
        end
    end
    
    return aj[1]
end
