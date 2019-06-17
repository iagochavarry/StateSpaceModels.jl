export simulate

"""
    simulate(ss::StateSpace, N::Int, S::Int)

Simulate S future scenarios up to N steps ahead. Returns a p x N x S matrix where the dimensions represent, respectively,
the number of series in the model, the number of steps ahead, and the number of scenarios.
"""
function simulate(ss::StateSpace, N::Int, S::Int)

    # Load estimated covariance matrices
    H = ss.covariance.H
    Q = ss.covariance.Q

    # Load system matrices
    n, p, m, r = size(ss.model)
    Z, T, R = ztr(ss.model)
    Z = Z[:, :, 1]

    # Load a, P, and F at last in-sample instant
    a0 = ss.smoother.alpha[end, :]
    P0 = ss.filter.P[:, :, end]
    F0 = ss.filter.F[:, :, end]
    
    # State and variance forecasts
    a = Matrix{Float64}(undef, N, m)
    P = Array{Float64, 3}(undef, m, m, N)
    F = Array{Float64, 3}(undef, p, p, N)

    # Probability distribution
    dist = Vector{Distribution}(undef, N)

    # Initialization
    a[1, :]    = T*a0
    P[:, :, 1]    = T*P0*T' + R*Q*R'
    F[:, :, 1]    = Z*P[:, :, 1]*Z' + H
    dist[1] = MvNormal(vec(Z*a[1, :]), Symmetric(F[:, :, 1]))
    sim = Array{Float64, 3}(undef, ss.model.dim.p, N, S)

    for t = 2:N
        a[t, :] = T*a[t-1, :]
        P[:, :, t] = T*P[:, :, t-1]*T' + R*Q*R'
        F[:, :, t] = Z*P[:, :, t]*Z' + H
        dist[t] = MvNormal(vec(Z*a[t, :]), Symmetric(F[:, :, t]))
    end

    for t = 1:N
        sim[:, t, :] = rand(dist[t], S)
    end

    return sim
end