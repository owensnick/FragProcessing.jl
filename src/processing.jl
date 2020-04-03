

gencsl(tbl; dw=1000, cf=:chrom, sf=:max_motif_start, ef=:max_motif_stop, strf=:max_motif_strand) = (tbl[cf], [div(a + b, 2) .+ (-dw:dw) for (a, b) in zip(tbl[sf], tbl[ef])], tbl[strf]);




@inline inc_meta_frag!(V, mpl, mpr, w) =  V[intersect(mpl:mpr, 1:size(V, 1))] .+= w

@inline inc_meta_cut!(V, mpl, mpr, w)  =  (inc_meta_index!(V, mpl, w) ; inc_meta_index!(V, mpr, w))
@inline inc_meta_mid!(V, mpl, mpr, w)  =  inc_meta_index!(V, div(mpl + mpr, 2), w)
@inline inc_meta_index!(V, i::Int, w)  =  (1 ≤ i ≤ size(V, 1)) && (V[i] += w)


@inline inc_heat_frag!(V, k, mpl, mpr, w) =  V[intersect(mpl:mpr, 1:size(V, 1)), k] .+= w
@inline inc_heat_cut!(V, k, mpl, mpr, w)  =  (inc_heat_index!(V, k, mpl, w) ; inc_heat_index!(V, k, mpr, w))
@inline inc_heat_left_cut!(V, k, mpl, mpr, w)  =  inc_heat_index!(V, k, mpl, w)
@inline inc_heat_mid!(V, k, mpl, mpr, w)  =  inc_heat_index!(V, k, div(mpl + mpr, 2), w)
@inline @inbounds inc_heat_index!(V, k, i::Int, w)  =  (1 ≤ i ≤ size(V, 1)) && (V[i, k] += w)


inc_meta_mid_width(bw) = (V, mpl, mpr, w) -> begin
    mid = div(mpl + mpr, 2)
    inc_meta_frag!(V, mid - bw, mid + bw, w)
end

inc_heat_mid_width(bw) = (V, k, mpl, mpr, w) -> begin
    mid = div(mpl + mpr, 2)
    bc = intersect((mid - bw):(mid + bw), 1:size(V, 1))
    V[bc, k] .+= w
end

inc_heat_atac_cut_width(bw) = (V, k, mpl, mpr, w, ao=4) -> begin
    inc_heat_frag!(V, k, mpl + ao - bw, mpl + ao + bw, w)
    inc_heat_frag!(V, k, mpr - ao - bw, mpr - ao + bw, w)
end

inc_meta_atac_cut_width(bw) = (V, mpl, mpr, w, ao=4) -> begin
    inc_meta_frag!(V, mpl + ao - bw, mpl + ao + bw, w)
    inc_meta_frag!(V, mpr - ao - bw, mpr - ao + bw, w)
end

@inline inc_heat_atac_cut!(V, k, mpl, mpr, w, ao=4)   =  (inc_heat_index!(V, k, mpl + ao, w) ; inc_heat_index!(V, k, mpr - ao, w))
@inline inc_heat_atac_frag!(V, k, mpl, mpr, w, ao=4)  =  inc_heat_frag!(V, k, mpl+ao, mpr-ao, w);
@inline inc_meta_atac_frag!(V, mpl, mpr, w, ao=4)     =  inc_meta_frag!(V, mpl+ao, mpr-ao, w);

function libsize(FM, dataentry, fn=identity, minfrag=0, maxfrag=1000000)

    σ = 0.0
    for i = 1:size(FM.FM, 2)
        fp = FM.FM[2, i] - FM.FM[1, i] + 1
        σ += ifelse(minfrag ≤ fp ≤ maxfrag, fn(FM.FM[dataentry, i]), 0.0)
    end
    σ
end



function libsize(FM, dataentry, fn=identity)
    σ = 0
    @inbounds for i = 1:size(FM.FM, 2)
        σ += fn(FM.FM[dataentry, i])
    end
    σ
end




function frag_heatmap_strand_norm(chroms, locations, strands, FM::Vector{FragMatrixSingle{T}}, inc_fun=inc_heat_mid!, fraglength=120, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    H = zeros(Float64, fw, length(chroms))
    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = zeros(Float64, fw, length(chroms))
        frag_heatmap_strand!(TH, p, chroms, locations, strands, F, inc_fun, fraglength, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
	    elseif norm_type == :frag_region
	        H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H += TH
        else
            H += TH
        end
    end
    if norm_type == :total
        rmul!(H, norm_scale/sum(F.totalecfrags[data_entry - 3] for F in FM))
    else
        rmul!(H, 1/length(FM))
    end
    H
end

function frag_heatmap_strand_norm(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    H = zeros(Float64, fw, length(chroms))
    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = zeros(Float64, fw, length(chroms))
        frag_heatmap_strand!(TH, p, chroms, locations, strands, F, inc_fun, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
        elseif norm_type == :frag
            σ = libsize(F, data_entry, fn, minfragsize, maxfragsize)
            H .+= norm_scale.*TH./σ
	elseif norm_type == :frag_region
	    H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H += TH
        else
            H += TH
        end
    end
    if norm_type == :total
        rmul!(H, 1e+6/sum(F.totalecfrags[data_entry - 3] for F in FM))
    else
        rmul!(H, 1/length(FM))
    end
    H
end


function frag_heatmap_strand(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])
    H = zeros(Float64, fw, length(chroms))
    p = Progress(length(chroms)*length(FM))
    for F in FM
        frag_heatmap_strand!(H, p, chroms, locations, strands, F, inc_fun, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
    end
    H
end

function frag_heatmap_strand(chroms, locations, strands, FM::Vector{FragMatrixSingle{T}}, inc_fun=inc_heat_mid!, fraglength=120, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])
    H = zeros(Float64, fw, length(chroms))
    p = Progress(length(chroms)*length(FM))
    for F in FM
        frag_heatmap_strand!(H, p, chroms, locations, strands, F, inc_fun, fraglength, data_entry, fn, show_progress=show_progress)
    end
    H
end

function frag_heatmap_strand!(H, p, chroms, locations, strands, FM::FragMatrixPair{T}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])

    for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, l, FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            fp = V[2, i] - V[1, i] + 1
            ((fp < minfragsize) || (fp > maxfragsize)) && continue
            if s == "+"
                mpl = V[1, i] - l.start + 1
                mpr = V[2, i] - l.start + 1
                strand = 1
            else
                mpl = l.stop - V[2, i] + 1
                mpr = l.stop - V[1, i] + 1
                strand = -1
            end

            ### flip missing?

            w = fn(V[data_entry, i])
            inc_fun(H, k, mpl, mpr,  w)
        end
    end
    H
end


function frag_heatmap_strand!(H, p, chroms, locations, strands, FM::FragMatrixSingle{T}, inc_fun=inc_heat_mid!, fraglength=120, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])

    for (k, (c, l, str)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, (l.start - 2*fraglength):(l.stop + 2*fraglength), FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            strand, _ = GenomeFragments.get_strand_chrom_enc(V[3, i])

            if strand == STRAND_POS
                s = V[1, i] - l.start + 1
                e = s + fraglength - 1
            else
                e = V[2, i]  - l.start + 1
                s = e - fraglength + 1
            end
            if str == "-"
                s, e = fw - e + 1, fw - s + 1
            end
            w = fn(V[data_entry, i])
            inc_fun(H, k, s, e,  w)
        end
    end
    H
end


function average_heatmap(H::Vector{T}, δ) where {T}
    n = length(H)

    w = cld(n, δ)
    AH = zeros(Float64, w)
    te = zeros(Int, w)
    for i = 1:n
        wi = cld(i, δ)
        te[wi] += 1
        AH[wi] += H[i]
    end
    if te[end] < 0.75δ
        AH[end] += AH[end-1]
        te[end] += te[end-1]
    end
    AH./te
end


average_heatmap(H, δh, δv) = average_heatmap_vert(average_heatmap(H, δh), δv)
max_heatmap(H, δh, δv) = max_heatmap_vert(max_heatmap(H, δh), δv)

function average_heatmap_vert(H, δ=2)

    n, m = size(H)

    w = cld(n, δ)
    AH = zeros(Float64, w, m)
    te = zeros(Int, w)
    for i = 1:m
        for j = 1:n
            wj = cld(j, δ)
            (i == 1) && (te[wj] += 1)
            AH[wj, i] += H[j, i]
        end
    end

    if te[end] < 0.75δ
        for i = 1:m
            AH[end, i] += AH[end-1, i]
        end
        te[end] += te[end-1]
    end

    AH./te
end

function average_heatmap(H, δ=2)

    n, m = size(H)
    w = cld(m, δ)

    AH = zeros(Float64, n, w)
    te = zeros(Int, w)
    for i = 1:m
        wi = cld(i, δ)
        te[wi] += 1
        for j = 1:n
            AH[j, wi] += H[j, i]
        end
    end
    if te[end] < 0.75δ
        for j = 1:n
            AH[j, end] += AH[j, end-1]
        end
        te[end] += te[end-1]
    end
    AH./te'
end


function max_heatmap_vert(H, δ=2)

    n, m = size(H)

    w = cld(n, δ)
    AH = zeros(Float64, w, m)
    te = zeros(Int, w)
    for i = 1:m
        for j = 1:n
            wj = cld(j, δ)
            (i == 1) && (te[wj] += 1)
            AH[wj, i] = max(AH[wj, i], H[j, i])
        end
    end

    if te[end] < 0.75δ
        for i = 1:m

            AH[end, i] = max(AH[end, i], AH[end-1, i])
        end
        te[end] += te[end-1]
    end
    AH
end


function max_heatmap(H, δ=2)

    n, m = size(H)
    w = cld(m, δ)

    AH = zeros(Float64, n, w)
    te = zeros(Int, w)
    for i = 1:m
        wi = cld(i, δ)
        te[wi] += 1
        for j = 1:n
            AH[j, wi] = max(AH[j, wi], H[j, i])
        end
    end
    if te[end] < 0.75δ
        for j = 1:n
            AH[j, end] = max(AH[j, end], AH[j, end-1])
        end
        te[end] += te[end-1]
    end
    # AH./te'
    AH
end



function frag_vplot_strand_norm(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    fr = length(minfragsize:maxfragsize)
    H = zeros(Float64, fw, fr)
    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = zeros(Float64, fw, fr)
        frag_vplot_strand!(TH, p, chroms, locations, strands, F, inc_fun, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
        elseif norm_type == :frag
            σ = libsize(F, data_entry, fn, minfragsize, maxfragsize)
            H .+= norm_scale.*TH./σ
	elseif norm_type == :frag_region
	    H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H += TH
        else
            H += TH
        end
    end
    if norm_type == :total
        rmul!(H, 1e+6/sum(F.totalecfrags[data_entry - 3] for F in FM))
    else
        rmul!(H, 1/length(FM))
    end
    H
end


function frag_vplot_strand(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])
    fr = length(minfragsize:maxfragsize)
    H = zeros(Float64, fw, fr)
    p = Progress(length(chroms)*length(FM))
    for F in FM
        frag_vplot_strand!(H, p, chroms, locations, strands, F, inc_fun, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
    end
    H
end

function frag_vplot_strand!(H, p, chroms, locations, strands, FM::FragMatrixPair{T}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])

    for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, l, FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            fp = V[2, i] - V[1, i] + 1
            ((fp < minfragsize) || (fp > maxfragsize)) && continue
            if s == "+"
                mpl = V[1, i] - l.start + 1
                mpr = V[2, i] - l.start + 1
                strand = 1
            else
                mpl = l.stop - V[2, i] + 1
                mpr = l.stop - V[1, i] + 1
                strand = -1
            end

            ### flip missing?

            w = fn(V[data_entry, i])
            inc_fun(H, fp, mpl, mpr,  w)
        end
    end
    H
end

function average_heatmap(H::Vector{T}, δ) where {T}
    n = length(H)

    w = cld(n, δ)
    AH = zeros(Float64, w)
    te = zeros(Int, w)
    for i = 1:n
        wi = cld(i, δ)
        te[wi] += 1
        AH[wi] += H[i]
    end
    if te[end] < 0.75δ
        AH[end] += AH[end-1]
        te[end] += te[end-1]
    end
    AH./te
end


#############################




function frag_meta_strand_norm(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    H = zeros(Float64, fw)
    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = zeros(Float64, fw)
        frag_meta_strand!(TH, p, chroms, locations, strands, F, inc_fun, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
        elseif norm_type == :frag
            σ = libsize(F, data_entry, fn, minfragsize, maxfragsize)
            H .+= norm_scale.*TH./σ
	elseif norm_type == :frag_region
	    H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H += TH
        else
            H += TH
        end
    end
    if norm_type == :total
        rmul!(H, 1e+6/sum(F.totalecfrags[data_entry - 3] for F in FM))
    else
        rmul!(H, 1/length(FM))
    end

    rmul!(H, 1/length(chroms))
    H
end


function frag_meta_strand(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])
    H = zeros(Float64, fw)
    p = Progress(length(chroms)*length(FM))
    for F in FM
        frag_meta_strand!(H, p, chroms, locations, strands, F, inc_fun, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
    end
    rmul!(H, 1/length(chroms))
    H
end


function frag_meta_strand!(H, p, chroms, locations, strands, FM::FragMatrixPair{T}, inc_fun=inc_heat_mid!, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])

    for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, l, FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            fp = V[2, i] - V[1, i] + 1
            ((fp < minfragsize) || (fp > maxfragsize)) && continue
            if s == "+"
                mpl = V[1, i] - l.start + 1
                mpr = V[2, i] - l.start + 1
                strand = 1
            else
                mpl = l.stop - V[2, i] + 1
                mpr = l.stop - V[1, i] + 1
                strand = -1
            end

            w = fn(V[data_entry, i])
            inc_fun(H, mpl, mpr,  w)
        end
    end
    H
end


function fragregion(chrom, location, strand, FM::FragMatrixSingle{T}, inc_fun=inc_meta_frag!, fraglength=120, dataentry=4, fn=identity; pos=true, neg=true) where {T}
    fw = length(location)
    P = zeros(fw)
    V = get_frags(chrom, location, FM)

    for i = 1:size(V, 2)
        rstrand, _ = GenomeFragments.get_strand_chrom_enc(V[3, i])
        !pos && (rstrand == STRAND_POS) && continue
        !neg && (rstrand == STRAND_NEG) && continue
        if rstrand == STRAND_POS
            s = V[1, i] - location.start + 1
            e = s + fraglength - 1
        else
            e = V[2, i]  - location.start + 1
            s = e - fraglength + 1
        end
        if strand == "-"
            s, e = fw - e + 1, fw - s + 1
        end
        w = fn(V[dataentry, i])
        inc_fun(P, s, e,  w)
    end
    P
end


@inline function fragregion(chrom, location, FM::FragMatrixPair{T}, inc_fun=inc_meta_frag!, minfrag=0, maxfrag=1000, dataentry=4, fn=identity) where {T}
    fw = length(location)
    P = zeros(fw)
    V = get_frags(chrom, location, FM)

    for col in eachcol(V)
        fp = col[2] - col[1] + 1
        s = col[1] - location.start + 1
        e = col[2] - location.start + 1
        w = fn(col[dataentry])
        (minfrag ≤ fp ≤ maxfrag) && inc_fun(P, s, e, w)
    end

    P
end

@inline function fragregion(chrom, location, FMS::Vector{FragMatrixPair{T}}, inc_fun=inc_meta_frag!, minfrag=0, maxfrag=1000, dataentry=4, fn=identity) where {T}
    fw = length(location)
    P = zeros(fw)

    for FM in FMS
        V = get_frags(chrom, location, FM)

        for col in eachcol(V)
            fp = col[2] - col[1] + 1
            s = col[1] - location.start + 1
            e = col[2] - location.start + 1
            w = fn(col[dataentry])
            (minfrag ≤ fp ≤ maxfrag) && inc_fun(P, s, e, w)
        end

        P
    end
    P
end



function count_cut_sites(chroms, locations, strands, FM, minfragsize=0, maxfragsize=1000, data_entry=4, fn=identity) where T

    cutsites = zeros(Int, length(chroms))
    for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))

        V = get_frags(c, l, FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            fp = V[2, i] - V[1, i] + 1
            ((fp < minfragsize) || (fp > maxfragsize)) && continue

            w = fn(V[data_entry, i])
            cutsites[k] += ifelse(V[1, i] ∈ l, w, 0)
            cutsites[k] += ifelse(V[2, i] ∈ l, w, 0)
        end
    end
    cutsites


end


function readsintersecting(chroms, locations, FM, dataentry=4, fn=identity)
    total = zeros(Int, length(chroms))
    for (k, (c, l)) in enumerate(zip(chroms, locations))

        v = get_frags(c, l, FM)
        for i = 1:size(v, 2)
            total[k] += fn(v[dataentry, i])
        end
    end
    total
end

function peakheights(chroms, locations, FM::Vector{FragMatrixPair{T}} ; minfragsize=0, maxfragsize=1000, dataentry=4, fn=one, scale=1e+6) where {T}

    PH = zeros(length(chroms))
    σ = 0
    σ = sum(libsize(F, dataentry, fn) for F in FM)

    for (k, (c, l)) in enumerate(zip(chroms, locations))
        PH[k] = maximum(mapreduce(F -> fragregion(c, l, F, inc_meta_frag!, minfragsize, maxfragsize, dataentry, fn), +, FM))
    end

    NH = scale.*PH./σ
    NH
end


function peakheightsloc(chroms, locations, FM ; minfragsize=0, maxfragsize=1000, dataentry=4, fn=one, scale=1e+6)

    PH = zeros(length(chroms))
    PS = zeros(Int, length(chroms))

    for (k, (c, l)) in enumerate(zip(chroms, locations))
        PH[k], PS[k] = findmax(fragregion(c, l, FM, inc_meta_frag!, minfragsize, maxfragsize, dataentry, fn))
    end

    PH, PS
end

function frag_heat_cut_pair_strand_norm(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    # fr = length(minfragsize:maxfragsize)
    H = zeros(Float64, fw, fw)
    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = zeros(Float64, fw, fw)
        frag_heat_cut_pair_strand!(TH, p, chroms, locations, strands, F, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
        elseif norm_type == :frag
            σ = libsize(F, data_entry, fn, minfragsize, maxfragsize)
            H .+= norm_scale.*TH./σ
	elseif norm_type == :frag_region
	    H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H += TH
        else
            H += TH
        end
    end
    if norm_type == :total
        rmul!(H, 1e+6/sum(F.totalecfrags[data_entry - 3] for F in FM))
    else
        rmul!(H, 1/length(FM))
    end
    H
end


function frag_heat_cut_pair_strand!(H, p, chroms, locations, strands, FM::FragMatrixPair{T}, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])

    for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, l, FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            fp = V[2, i] - V[1, i] + 1
            ((fp < minfragsize) || (fp > maxfragsize)) && continue
            if s == "+"
                mpl = V[1, i] - l.start + 1
                mpr = V[2, i] - l.start + 1
                strand = 1
            else
                mpl = l.stop - V[2, i] + 1
                mpr = l.stop - V[1, i] + 1
                strand = -1
            end

            w = fn(V[data_entry, i])
            if (1 <= mpl <= fw) && (1 <= mpr <= fw)
                H[mpl, mpr] += w
            end
        end
    end
    H
end



function frag_heat_cut_pair_strand_norm(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    # fr = length(minfragsize:maxfragsize)
    H = zeros(Float64, fw, fw)
    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = zeros(Float64, fw, fw)
        frag_heat_cut_pair_strand!(TH, p, chroms, locations, strands, F, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
        elseif norm_type == :frag
            σ = libsize(F, data_entry, fn, minfragsize, maxfragsize)
            H .+= norm_scale.*TH./σ
	elseif norm_type == :frag_region
	    H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H += TH
        else
            H += TH
        end
    end
    if norm_type == :total
        rmul!(H, 1e+6/sum(F.totalecfrags[data_entry - 3] for F in FM))
    else
        rmul!(H, 1/length(FM))
    end
    H
end


function frag_tensor_cut_pair_strand_norm(chroms, locations, strands, FM::Vector{FragMatrixPair{T}}, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; norm_scale=1e+6, norm_type=:mpm, show_progress=true) where T

    fw = length(locations[1])
    fr = length(chroms)
    H = [spzeros(Float64, fw, fw) for i = 1:fr]

    p = ProgressMeter.Progress(length(chroms)*length(FM))
    for F in FM
        TH = [spzeros(Float64, fw, fw) for i = 1:fr]
        frag_tensor_cut_pair_strand!(TH, p, chroms, locations, strands, F, minfragsize, maxfragsize, data_entry, fn, show_progress=show_progress)
        if norm_type == :mpm
            H .+= norm_scale.*TH./F.totalecfrags[data_entry - 3]
        elseif norm_type == :mpm_f
            t = sum(fn, view(F.FM, data_entry, :))
            H .+= norm_scale.*TH./t
        elseif norm_type == :frag
            σ = libsize(F, data_entry, fn, minfragsize, maxfragsize)
            H .+= norm_scale.*TH./σ
	    elseif norm_type == :frag_region
	        H .+= norm_scale.*TH/sum(TH)
        elseif norm_type == :total
            H .+= TH
        else
            H .+= TH
        end
    end
    if norm_type == :total
        #rmul!(H,)
        f = 1e+6/sum(F.totalecfrags[data_entry - 3] for F in FM)
        H .*= f
    else
        H ./= length(FM)
    end
    H
end

function frag_tensor_cut_pair_strand!(H, p, chroms, locations, strands, FM::FragMatrixPair{T}, minfragsize=0, maxfragsize=500, data_entry=4, fn=identity; show_progress=true) where T

    fw = length(locations[1])

    for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))
        show_progress && next!(p)
        V = get_frags(c, l, FM)
        isempty(V) && continue
        for i = 1:size(V, 2)
            fp = V[2, i] - V[1, i] + 1
            ((fp < minfragsize) || (fp > maxfragsize)) && continue
            if s == "+"
                mpl = V[1, i] - l.start + 1
                mpr = V[2, i] - l.start + 1
                strand = 1
            else
                mpl = l.stop - V[2, i] + 1
                mpr = l.stop - V[1, i] + 1
                strand = -1
            end

            w = fn(V[data_entry, i])
            if (1 <= mpl <= fw) && (1 <= mpr <= fw)
                H[k][mpl, mpr] += w
            end
        end
    end
    H
end
