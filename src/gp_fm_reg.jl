
using PyCall, Statistics, LinearAlgebra

gpy = pyimport("GPy")
gpy_kern = pyimport("GPy.kern")
gpy_models = pyimport("GPy.models")


data_trans(X, α=1e+6, β=1) = log(α*X + β)
data_inv_trans(X, α=1e+6, β=1) = (exp(X) - β)/α
spectraldensity(s, l) = l*sqrt(2*π)*exp(-2*(π*s*l).^2)

get_train_ind(X, train_region) =  abs.(X) .<= train_region
function get_train_ind(X, train_region::T) where T <: AbstractVector
    first(train_region) .<= abs.(X) .<= last(train_region)
end


function specdensity(nxp, X, nr=6, ; trans=data_trans, inv_trans=data_inv_trans, postnorm=identity, train_region=100:600)
    any(isnan, X) && return (0.0, 0.0, 0.0, 0.0), zeros(size(X))
    gp_p, train, pred, sfs  = gpy_reg_train_region_post_sample(nxp, postnorm(X), train_region, trans, inv_trans, post_samples=2)
    SF = vec(first(pred))
    return gp_p, SF
end


function gpy_reg_train_region_post_sample(SX, ρ, train_region=Inf, trans=data_trans, inv_trans=data_inv_trans ;  post_samples=10, full_cov=false)

    p, train, pred, Λ  = gpy_reg_sparse_train_post_sample(SX, trans.(ρ), train_region, 10, 10, 1, 1:100, post_samples, full_cov)

    X_train, μ_train, σ_train, ql_train, qu_train = train
    sf_train = inv_trans.(μ_train)
    cil_train = inv_trans.(ql_train)
    ciu_train = inv_trans.(qu_train)


    μ_pred, σ_pred, ql_pred, qu_pred = pred
    sf_pred  = inv_trans.(μ_pred)
    cil_pred = inv_trans.(ql_pred)
    ciu_pred = inv_trans.(qu_pred)
    sfs_pred = inv_trans.(Λ)


    p, (X_train, sf_train, cil_train, ciu_train), (sf_pred, cil_pred, ciu_pred, σ_pred), sfs_pred
end


function gpy_reg_sparse_train_post_sample(X, y, train_region, zstep_train=10, zstep_pred=10, num_restarts=0, var_region=1:100, post_samples=10, full_cov=false)

    𝕍y = var(y[var_region])
    𝕍f = var(y)
    initℓ = 140


    train_ind = get_train_ind(X, train_region)
    Xtrain = X[train_ind, 1:1]
    ytrain = y[train_ind, 1:1]
    Ztrain = unique(Xtrain)[1:zstep_train:end, 1:1]

    kern = gpy_kern.RBF(input_dim=1, variance=𝕍f, lengthscale=initℓ)
    gpr = gpy_models.SparseGPRegression(Xtrain, ytrain, kern, Z=Ztrain)
    gpr.Gaussian_noise.variance = 𝕍y

    gpr.optimize()
    (num_restarts > 0) && gpr.optimize_restarts(num_restarts=num_restarts)

    σn² = gpr.Gaussian_noise.variance[1]
    ℓ   = gpr.kern.lengthscale[1]
    σf² = gpr.kern.variance[1]


    μtrain, σtrain = gpr.predict(Xtrain)
    ciltrain, ciutrain = gpr.predict_quantiles(Xtrain)


    #### Now wider region
    X = X[:, 1:1]
    y = y[:, 1:1]
    Z = X[1:zstep_pred:end, 1:1]
    kern = gpy_kern.RBF(input_dim=1, variance=σf², lengthscale=ℓ)
    gpr = gpy_models.SparseGPRegression(X, y, kern, Z=Z)
    gpr.Gaussian_noise.variance = σn²
    lml = gpr.log_likelihood()[1]

    μ, σ = gpr.predict(X, full_cov=full_cov)
    cil, ciu = gpr.predict_quantiles(X)
    Λ = gpr.posterior_samples_f(X, post_samples)


    (σf², ℓ, σn², lml), (vec(Xtrain), μtrain, σtrain, ciltrain, ciutrain), (μ, σ, cil, ciu), Λ
end

function gpy_reg_sparse_fix_ell(X, y, train_region; ell=Float64[], ldelta=10, lstep=2, ldef=range(60, 500, length=20),  zstep_train=10, zstep_pred=10, num_restarts=0, var_region=1:100, post_samples=10, full_cov=false)


    𝕍y = var(y[var_region])
    𝕍f = var(y)
    initℓ = 140


    train_ind = get_train_ind(X, train_region)
    Xtrain = X[train_ind, 1:1]
    ytrain = y[train_ind, 1:1]
    Ztrain = unique(Xtrain)[1:zstep_train:end, 1:1]

    kern = gpy_kern.RBF(input_dim=1, variance=𝕍f, lengthscale=initℓ)
    gpr = gpy_models.SparseGPRegression(Xtrain, ytrain, kern, Z=Ztrain)
    gpr.Gaussian_noise.variance = 𝕍y

    gpr.optimize()
    (num_restarts > 0) && gpr.optimize_restarts(num_restarts=num_restarts)

    σn² = gpr.Gaussian_noise.variance[1]
    ℓ   = gpr.kern.lengthscale[1]
    σf² = gpr.kern.variance[1]


    μtrain, σtrain = gpr.predict(Xtrain)
    ciltrain, ciutrain = gpr.:predict_quantiles(Xtrain)


    ### Sample ell
    if isempty(ell)
        ell = sort([range(ℓ - sqrt(ℓ), ℓ + sqrt(ℓ), length=20) ; ldef])
    end

    lml_ell = zeros(length(ell))
    for (i, l) in enumerate(ell)
        py"$(gpr).kern.lengthscale.fix($(l))"
        gpr.optimize()
        lml_ell[i] = gpr.log_likelihood()[1]
    end

    #### Now wider region
    X = X[:, 1:1]
    y = y[:, 1:1]
    Z = X[1:zstep_pred:end, 1:1]
    kern = gpy_kern.RBF(input_dim=1, variance=σf², lengthscale=ℓ)
    gpr = gpy_models.SparseGPRegression(X, y, kern, Z=Z)
    gpr.Gaussian_noise.variance = σn²
    lml = gpr.log_likelihood()[1]

    μ, σ = gpr.predict(X, full_cov=full_cov)
    cil, ciu = gpr.predict_quantiles(X)



    (σf², ℓ, σn², lml), (vec(Xtrain), μtrain, σtrain, ciltrain, ciutrain), (μ, σ, cil, ciu), (ell, lml_ell)
end



function specdensity_ell(nxp, X, nr=6, ; ell = Float64[], trans=data_trans, inv_trans=data_inv_trans, postnorm=identity, train_region=100:600)

    any(isnan, X) && return (0.0, 0.0, 0.0, 0.0), zeros(size(X)), (Float64[], Float64[])

    p, train, pred, ell_stats  = gpy_reg_sparse_fix_ell(nxp, trans.(X), train_region, num_restarts=1, ell=ell)

    μ_pred, σ_pred, ql_pred, qu_pred = pred
    sf_pred  = inv_trans.(μ_pred)

    return p, sf_pred, ell_stats
end

function specdensity_gpp(nxp, X; trans=data_trans, inv_trans=data_inv_trans, postnorm=identity, train_region=100:600, p=1/180)

    any(isnan, X) && return (0.0, 0.0, 0.0, 0.0), zeros(size(X)), (Float64[], Float64[])

    σf², ℓ, σn², lml = gpy_reg_sparse_train(nxp, trans.(X), train_region, num_restarts=1)

    (sigf2=σf², ell=ℓ, sign2=σn², lml=lml, sp=spectraldensity(p, ℓ))
end




function gpy_reg_sparse_train(X, y, train_region; zstep_train=10, zstep_pred=10, num_restarts=0, var_region=1:100)

    𝕍y = var(y[var_region])
    𝕍f = var(y)
    initℓ = 140


    train_ind = get_train_ind(X, train_region)
    Xtrain = X[train_ind, 1:1]
    ytrain = y[train_ind, 1:1]
    Ztrain = unique(Xtrain)[1:zstep_train:end, 1:1]

    kern = gpy_kern.RBF(input_dim=1, variance=𝕍f, lengthscale=initℓ)
    gpr = gpy_models.SparseGPRegression(Xtrain, ytrain, kern, Z=Ztrain)
    gpr.Gaussian_noise.variance = 𝕍y

    gpr.optimize()
    (num_restarts > 0) && gpr.optimize_restarts(num_restarts=num_restarts)

    σn² = gpr.Gaussian_noise.variance[1]
    ℓ   = gpr.kern.lengthscale[1]
    σf² = gpr.kern.variance[1]
    lml = gpr.log_likelihood()[1]

    (σf², ℓ, σn², lml)
end
