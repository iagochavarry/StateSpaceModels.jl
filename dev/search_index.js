var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#StateSpaceModels.jl-Documentation-1",
    "page": "Home",
    "title": "StateSpaceModels.jl Documentation",
    "category": "section",
    "text": "StateSpaceModels.jl is a package for modeling, forecasting, and simulating time series in a state-space framework. Implementations were made based on the book \"Time Series Analysis by State Space Methods\" (2012) by James Durbin and Siem Jan Koopman. The notation of the variables in the code also follows the book."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "This package is registered in METADATA so you can Pkg.add it as follows:pkg> add StateSpaceModels"
},

{
    "location": "#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "Current features:Kalman filter and smoother\nSquare-root filter and smoother\nMaximum likelihood estimation\nForecasting\nMonte Carlo simulation\nMultivariate modeling\nUser-defined models (input any Z, T, and R)\nSeveral predefined models, including:\nBasic structural model (trend, slope, seasonal)\nStructural model with exogenous variables\nLinear trend model\nLocal level model\nCompletion of missing values\nDiagnostics for the residuals, including:\nJarque-Bera test\nLjung-Box test\nHomoscedasticity testPlanned features:Exact initialization of the Kalman filter\nEM algorithm for maximum likelihood estimation\nUnivariate treatment of multivariate models"
},

{
    "location": "#Citing-StateSpaceModels.jl-1",
    "page": "Home",
    "title": "Citing StateSpaceModels.jl",
    "category": "section",
    "text": "If you use StateSpaceModels.jl in your work, we kindly ask you to cite the following paper (pdf):@article{SaavedraBodinSouto2019,\ntitle={StateSpaceModels.jl: a Julia Package for Time-Series Analysis in a State-Space Framework},\nauthor={Raphael Saavedra and Guilherme Bodin and Mario Souto},\njournal={arXiv preprint arXiv:1908.01757},\nyear={2019}\n}"
},

{
    "location": "manual/#",
    "page": "Manual",
    "title": "Manual",
    "category": "page",
    "text": ""
},

{
    "location": "manual/#Manual-1",
    "page": "Manual",
    "title": "Manual",
    "category": "section",
    "text": ""
},

{
    "location": "manual/#Introduction-1",
    "page": "Manual",
    "title": "Introduction",
    "category": "section",
    "text": "In this package we consider the following state-space modelbegingather*\n    beginaligned\n        y_t = Z_t alpha_t + d_t + varepsilon_t quad quad quad t = 1 dots n \n        alpha_t+1 = T alpha_t + c_t + R eta_t\n    endaligned\nendgather*wherebeginbmatrix\n    varepsilon_t \n    eta_t \n    alpha_1\nendbmatrix\nsim\nNID\nbeginpmatrix\n    beginbmatrix\n        0 \n        0 \n        a_1\n    endbmatrix\n    \n    beginbmatrix\n        H  0  0\n        0  Q  0\n        0  0  P_1\n    endbmatrix\nendpmatrix"
},

{
    "location": "manual/#StateSpaceModels.StateSpaceDimensions",
    "page": "Manual",
    "title": "StateSpaceModels.StateSpaceDimensions",
    "category": "type",
    "text": "StateSpaceDimensions\n\nStateSpaceModel dimensions, following the notation of on the book \"Time Series Analysis by State Space Methods\" (2012) by J. Durbin and S. J. Koopman.\n\nn is the number of observations\np is the dimension of the observation vector y_t\nm is the dimension of the state vector alpha_t\nr is the dimension of the state covariance matrix Q\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.StateSpaceModel",
    "page": "Manual",
    "title": "StateSpaceModels.StateSpaceModel",
    "category": "type",
    "text": "StateSpaceModel{Typ}\n\nFollowing the notation of on the book \"Time Series Analysis by State Space Methods\" (2012) by J. Durbin and S. J. Koopman.\n\ny A n times p matrix containing observations\nZ A p times m times n matrix\nT A m times m matrix\nR A m times r matrix\nd A n times p matrix\nc A n times m matrix\nH A p times p matrix\nQ A r times r matrix\n\nThere are multiple constructors for the StateSpaceModel:\n\nStateSpaceModel(y::VecOrMat{Typ}, Z::Matrix{Typ}, T::Matrix{Typ}, R::Matrix{Typ}) where Typ <: Real StateSpaceModel(y::VecOrMat{Typ}, Z::Array{Typ, 3}, T::Matrix{Typ}, R::Matrix{Typ}) where Typ <: Real StateSpaceModel(y::VecOrMat{Typ}, Z::Matrix{Typ}, T::Matrix{Typ}, R::Matrix{Typ}, d::Matrix{Typ}, c::Matrix{Typ}, H::Matrix{Typ}, Q::Matrix{Typ}) where Typ <: Real StateSpaceModel(y::VecOrMat{Typ}, Z::Array{Typ, 3}, T::Matrix{Typ}, R::Matrix{Typ}, d::Matrix{Typ}, c::Matrix{Typ}, H::Matrix{Typ}, Q::Matrix{Typ}) where Typ <: Real\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.SmoothedState",
    "page": "Manual",
    "title": "StateSpaceModels.SmoothedState",
    "category": "type",
    "text": "SmoothedState\n\nFollowing the notation of on the book \"Time Series Analysis by State Space Methods\" (2012) by J. Durbin and S. J. Koopman.\n\nalpha Expected value of the smoothed state E(alpha_ty_1 dots  y_n)\nV Error covariance matrix of smoothed state Var(alpha_ty_1 dots  y_n)\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.FilterOutput",
    "page": "Manual",
    "title": "StateSpaceModels.FilterOutput",
    "category": "type",
    "text": "FilterOutput\n\nFollowing the notation of on the book \"Time Series Analysis by State Space Methods\" (2012) by J. Durbin and S. J. Koopman.\n\na Predictive state E(alpha_ty_t-1 dots  y_1)\natt Filtered state E(alpha_ty_t dots  y_1)\nv Prediction error v_t = y_t  Z_t a_t\nP Covariance matrix of predictive state P = Var(alpha_ty_t1 dots  y_1)\nPtt Covariance matrix of filtered state P = Var(alpha_ty_t dots  y_1)\nF Variance of prediction error Var(v_t)\nsteadystate Boolean to indicate if steady state was attained\ntsteady Instant when steady state was attained; in case it was not, tsteady = n+1\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.StateSpace",
    "page": "Manual",
    "title": "StateSpaceModels.StateSpace",
    "category": "type",
    "text": "StateSpace\n\nA state-space structure containing the model, filter output, smoother output, covariance matrices, filter type and optimization method.\n\n\n\n\n\n"
},

{
    "location": "manual/#Data-structures-1",
    "page": "Manual",
    "title": "Data structures",
    "category": "section",
    "text": "StateSpaceDimensions\nStateSpaceModel\nSmoothedState\nFilterOutput\nStateSpace"
},

{
    "location": "manual/#StateSpaceModels.local_level",
    "page": "Manual",
    "title": "StateSpaceModels.local_level",
    "category": "function",
    "text": "local_level(y::VecOrMat{Typ}) where Typ <: Real\n\nBuild state-space system for a local level model with observations y.\n\nIf y is provided as an Array{Typ, 1} it will be converted to an Array{Typ, 2} inside the StateSpaceModel.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.linear_trend",
    "page": "Manual",
    "title": "StateSpaceModels.linear_trend",
    "category": "function",
    "text": "linear_trend(y::VecOrMat{Typ}) where Typ <: Real\n\nBuild state-space system for a linear trend model with observations y.\n\nIf y is provided as an Array{Typ, 1} it will be converted to an Array{Typ, 2} inside the StateSpaceModel.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.structural",
    "page": "Manual",
    "title": "StateSpaceModels.structural",
    "category": "function",
    "text": "structural(y::VecOrMat{Typ}, s::Int; X::VecOrMat{Typ} = Matrix{Typ}(undef, 0, 0)) where Typ <: Real\n\nBuild state-space system for a given structural model with observations y, seasonality s, and, optionally, exogenous variables X.\n\nIf y is provided as an Array{Typ, 1} it will be converted to an Array{Typ, 2} inside the StateSpaceModel. The same will happen to X,  if an Array{Typ, 1} it will be converted to an Array{Typ, 2} inside the StateSpaceModel.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.regression",
    "page": "Manual",
    "title": "StateSpaceModels.regression",
    "category": "function",
    "text": "regression(y::VecOrMat{Typ}, X::VecOrMat{Typ}) where Typ <: Real\n\nBuild state-space system for estimating a regression model y_t = X_tbeta_t + varpsilon_t. Once the model is estimated the user can recover the parameter hat beta by the querying the smoothed states of the model.\n\n\n\n\n\n"
},

{
    "location": "manual/#Predefined-models-1",
    "page": "Manual",
    "title": "Predefined models",
    "category": "section",
    "text": "The package provides constructors for some classic state-space models.The local level model is defined bybegingather*\n    beginaligned\n        y_t =  mu_t  + varepsilon_t quad varepsilon_t sim mathcalN(0 sigma^2_varepsilon)\n        mu_t+1 = mu_t + eta_t quad eta_t sim mathcalN(0 sigma^2_eta)\n    endaligned\nendgather*local_levelThe linear trend model is defined bybegingather*\n    beginaligned\n        y_t =  mu_t  + varepsilon_t quad varepsilon_t sim mathcalN(0 sigma^2_varepsilon)\n        mu_t+1 = mu_t + nu_t + xi_t quad xi_t sim mathcalN(0 sigma^2_xi)\n        nu_t+1 = nu_t + zeta_t quad zeta_t sim mathcalN(0 sigma^2_zeta)\n    endaligned\nendgather*linear_trendThe structural model is defined bybegingather*\n    beginaligned\n        y_t =  mu_t + gamma_t + varepsilon_t quad varepsilon_t sim mathcalN(0 sigma^2_varepsilon)\n        mu_t+1 = mu_t + nu_t + xi_t quad xi_t sim mathcalN(0 sigma^2_xi)\n        nu_t+1 = nu_t + zeta_t quad zeta_t sim mathcalN(0 sigma^2_zeta)\n        gamma_t+1 = sum_j=1^s-1 gamma_t+1-j + omega_t quad  omega_t sim mathcalN(0 sigma^2_omega)\n    endaligned\nendgather*structuralA regression can also be defined in terms of state-space model. We consider the simple model y_t = X_tbeta_t + varpsilon_t where  varepsilon_t sim mathcalN(0 H_t). This can be written in the state-space form by doing Z_t = X_t T_t = I R = 0 and Q = 0regression"
},

{
    "location": "manual/#Custom-models-1",
    "page": "Manual",
    "title": "Custom models",
    "category": "section",
    "text": "You can also build your own custom model by passing all or some of the matrices that compose the state-space model.Currently there is a list of constructors that allow you to define a new model. When building the StateSpaceModel the parameters to be estimated will be indicated with a NaN. StateSpaceModel(y::Matrix{Typ}, Z::Array{Typ, 3}, T::Matrix{Typ}, R::Matrix{Typ}, H::Matrix{Typ}, Q::Matrix{Typ}) where Typ <: Real\n StateSpaceModel(y::Matrix{Typ}, Z::Matrix{Typ}, T::Matrix{Typ}, R::Matrix{Typ}, H::Matrix{Typ}, Q::Matrix{Typ}) where Typ <: Real\n StateSpaceModel(y::Matrix{Typ}, Z::Array{Typ, 3}, T::Matrix{Typ}, R::Matrix{Typ}) where Typ <: Real\n StateSpaceModel(y::Matrix{Typ}, Z::Matrix{Typ}, T::Matrix{Typ}, R::Matrix{Typ}) where Typ <: RealIn case H and Q are not provided, they are automatically filled with NaN and set to be estimated."
},

{
    "location": "manual/#StateSpaceModels.statespace",
    "page": "Manual",
    "title": "StateSpaceModels.statespace",
    "category": "function",
    "text": "statespace(model::StateSpaceModel; filter_type::DataType = KalmanFilter, optimization_method::AbstractOptimizationMethod = RandomSeedsLBFGS(), verbose::Int = 1)\n\nEstimate the pre-specified state-space model. The function will only estimate the entries that are declared NaN. If there are no NaNs in the model it will only perform the filter and smoother computations.\n\n\n\n\n\n"
},

{
    "location": "manual/#Estimation-1",
    "page": "Manual",
    "title": "Estimation",
    "category": "section",
    "text": "The model estimation is made using the function statespace(model; filter_type = KalmanFilter, optimization_method = RandomSeedsLBFGS(), verbose = 1). It receives as argument the pre-specified StateSpaceModel object model. Optionally, the user can define the Kalman filter variant to be used, the optimization method and the verbosity level.statespace"
},

{
    "location": "manual/#StateSpaceModels.forecast",
    "page": "Manual",
    "title": "StateSpaceModels.forecast",
    "category": "function",
    "text": "forecast(ss::StateSpace{Typ}, N::Int) where Typ\n\nObtain the minimum mean square error forecasts N steps ahead. Returns the forecasts and the predictive distributions  at each time period.\n\n\n\n\n\n"
},

{
    "location": "manual/#Forecasting-1",
    "page": "Manual",
    "title": "Forecasting",
    "category": "section",
    "text": "Forecasting is conducted with the function forecast. It receives as argument a StateSpace object and the number of steps ahead N.forecast"
},

{
    "location": "manual/#StateSpaceModels.simulate",
    "page": "Manual",
    "title": "StateSpaceModels.simulate",
    "category": "function",
    "text": "simulate(ss::StateSpace{Typ}, N::Int, S::Int) where Typ\n\nSimulate S future scenarios up to N steps ahead. Returns a p x N x S matrix where the dimensions represent, respectively, the number of series in the model, the number of steps ahead, and the number of scenarios.\n\n\n\n\n\n"
},

{
    "location": "manual/#Simulation-1",
    "page": "Manual",
    "title": "Simulation",
    "category": "section",
    "text": "Simulation is made using the function simulate. It receives as argument a StateSpace object, the number of steps ahead N and the number of scenarios to simulate S.simulate"
},

{
    "location": "manual/#Filters-1",
    "page": "Manual",
    "title": "Filters",
    "category": "section",
    "text": ""
},

{
    "location": "manual/#StateSpaceModels.kalman_filter",
    "page": "Manual",
    "title": "StateSpaceModels.kalman_filter",
    "category": "function",
    "text": "kalman_filter(model::StateSpaceModel{Typ}; tol::Typ = Typ(1e-5)) where Typ\n\nKalman filter with big Kappa initialization, i.e., initializing state variances as 1e6.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.smoother",
    "page": "Manual",
    "title": "StateSpaceModels.smoother",
    "category": "function",
    "text": "smoother(model::StateSpaceModel{Typ}, kfilter::KalmanFilter{Typ}) where Typ\n\nSmoother for state-space model.\n\n\n\n\n\n"
},

{
    "location": "manual/#Kalman-Filter-1",
    "page": "Manual",
    "title": "Kalman Filter",
    "category": "section",
    "text": "The implementation of the Kalman Filter follows the recursionbegingather*\n    beginaligned\n        v_t = y_t - Z_t a_t  F_t = Z_t P_t Z^top_t + H\n        a_tt = a_t - P_t Z^top_t F_t^-1 v_t quad quad quad P_tt = P_t -  P_t Z^top_t F_t^-1Z_t P_t\n        a_t+1 = T a_t - K_t v_t quad quad quad P_t+1 = T P_t(T - K_t Z_t)^top + R Q R\n    endaligned\nendgather*where K_t = T P_t Z^top_t F_t^-1. The terms a_t+1 and P_t+1 can be simplifed tobegingather*\n    beginaligned\n        a_t+1 = T a_ttquad quad quad P_t+1 = T P_tt T^top + R Q R\n    endaligned\nendgather*In case of missing observation the mean of inovations v_t and variance of inovations F_t become NaN and the recursion becomesbegingather*\n    beginaligned\n        a_tt = a_t quad quad quad P_tt = P_t\n        a_t+1 = T a_t quad quad quad P_t+1 = T P_t T^top + R Q R\n    endaligned\nendgather*StateSpaceModels.kalman_filterThe implementation of the smoother follows the recursionbegingather*\n    beginaligned\n        r_t-1 = Z^top_t F_t^-1 v_t + L^top_t r_t quad quad  N_t-1 = Z^top_t F_t^-1 Z_t + L^top_t N_t L_t\n        alpha_t = a_t + P_t r_t  V_t = P_t - P_t N_t-1 P_t  \n    endaligned\nendgather*where L_t = T - K_t Z_t.In case of missing observation then r_t-1 and N_t-1 becomebegingather*\n    beginaligned\n        r_t-1 = T^top r_t quad quad  N_t-1 = T^top_t N_t T_t\n    endaligned\nendgather*StateSpaceModels.smoother"
},

{
    "location": "manual/#StateSpaceModels.sqrt_kalman_filter",
    "page": "Manual",
    "title": "StateSpaceModels.sqrt_kalman_filter",
    "category": "function",
    "text": "sqrt_kalman_filter(model::StateSpaceModel{Typ}; tol::Typ = Typ(1e-5)) where Typ\n\nSquare-root Kalman filter with big Kappa initialization.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.sqrt_smoother",
    "page": "Manual",
    "title": "StateSpaceModels.sqrt_smoother",
    "category": "function",
    "text": "sqrt_smoother(model::StateSpaceModel, sqrt_filter::SquareRootFilter) where Typ\n\nSquare-root smoother for state space model.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.filtered_state",
    "page": "Manual",
    "title": "StateSpaceModels.filtered_state",
    "category": "function",
    "text": "filtered_state(model::StateSpaceModel{Typ}, sqrt_filter::SquareRootFilter{Typ}) where Typ\n\nObtain the filtered state estimates and their covariance matrices.\n\n\n\n\n\n"
},

{
    "location": "manual/#Square-Root-Kalman-Filter-1",
    "page": "Manual",
    "title": "Square Root Kalman Filter",
    "category": "section",
    "text": "StateSpaceModels.sqrt_kalman_filter\nStateSpaceModels.sqrt_smoother\nStateSpaceModels.filtered_state"
},

{
    "location": "manual/#StateSpaceModels.AbstractFilter",
    "page": "Manual",
    "title": "StateSpaceModels.AbstractFilter",
    "category": "type",
    "text": "AbstractFilter\n\nAbstract type used to implement an interface for generic filters.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.AbstractSmoother",
    "page": "Manual",
    "title": "StateSpaceModels.AbstractSmoother",
    "category": "type",
    "text": "AbstractSmoother\n\nAbstract type used to implement an interface for generic smoothers.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.kfas",
    "page": "Manual",
    "title": "StateSpaceModels.kfas",
    "category": "function",
    "text": "kfas(model::StateSpaceModel{T}, filter_type::DataType) where T\n\nPerform Kalman filter and smoother according to the chosen filter_type.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.statespace_likelihood",
    "page": "Manual",
    "title": "StateSpaceModels.statespace_likelihood",
    "category": "function",
    "text": "statespacelikelihood(psitilde::Vector{T}, model::StateSpaceModel, unknowns::Unknowns,                              validinsts::Vector{Int}, filter_type::DataType) where T\n\nCompute log-likelihood concerning hyperparameter vector psitilde (psi)\n\nEvaluate ell(psiy_n)= -fracnp2log2pi - frac12 sum_t=1^n log F_t -  frac12 sum_t=1^n v_t^top F_t^-1 v_t\n\n\n\n\n\n"
},

{
    "location": "manual/#Filter-interface-1",
    "page": "Manual",
    "title": "Filter interface",
    "category": "section",
    "text": "StateSpaceModels has an interface that allows users to define their own versions of the Kalman filter.StateSpaceModels.AbstractFilter\nStateSpaceModels.AbstractSmoother\nkfasEvery filter must provide its own version of the statespace_likelihood functionStateSpaceModels.statespace_likelihood"
},

{
    "location": "manual/#Optimization-methods-1",
    "page": "Manual",
    "title": "Optimization methods",
    "category": "section",
    "text": "StateSpaceModels has an interface that allows users to define their own optimization methods. It is easily integrated with Optim.jl."
},

{
    "location": "manual/#StateSpaceModels.LBFGS",
    "page": "Manual",
    "title": "StateSpaceModels.LBFGS",
    "category": "type",
    "text": "LBFGS(model::StateSpaceModel{T}, args...; kwargs...)\n\nIf an Int is provided the method will sample random seeds and use them as initial points for Optim LBFGS method. If a Vector{Vector{T}} is provided it will use them as initial points for Optim LBFGS method.\n\n\n\n\n\n"
},

{
    "location": "manual/#LBFGS-1",
    "page": "Manual",
    "title": "LBFGS",
    "category": "section",
    "text": "StateSpaceModels.LBFGS"
},

{
    "location": "manual/#StateSpaceModels.BFGS",
    "page": "Manual",
    "title": "StateSpaceModels.BFGS",
    "category": "type",
    "text": "BFGS(model::StateSpaceModel{T}, args...; kwargs...)\n\nIf an Int is provided the method will sample random seeds and use them as initial points for Optim BFGS method. If a Vector{Vector{T}} is provided it will use them as initial points for Optim BFGS method.\n\n\n\n\n\n"
},

{
    "location": "manual/#BFGS-1",
    "page": "Manual",
    "title": "BFGS",
    "category": "section",
    "text": "StateSpaceModels.BFGS"
},

{
    "location": "manual/#StateSpaceModels.AbstractOptimizationMethod",
    "page": "Manual",
    "title": "StateSpaceModels.AbstractOptimizationMethod",
    "category": "type",
    "text": "AbstractOptimizationMethod{T}\n\nAbstract type used to implement an interface for generic optimization methods.\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.estimate_statespace",
    "page": "Manual",
    "title": "StateSpaceModels.estimate_statespace",
    "category": "function",
    "text": "estimate_statespace(model::StateSpaceModel{T}, filter_type::DataType,\n                        optimization_method::AbstractOptimizationMethod; verbose::Int = 1) where T\n\nEstimate parameters of the StateSpaceModel according to its filter_type and optimization_method.\n\n\n\n\n\n"
},

{
    "location": "manual/#Optimization-methods-interface-1",
    "page": "Manual",
    "title": "Optimization methods interface",
    "category": "section",
    "text": "StateSpaceModels.AbstractOptimizationMethod\nStateSpaceModels.estimate_statespace"
},

{
    "location": "manual/#StateSpaceModels.diagnostics",
    "page": "Manual",
    "title": "StateSpaceModels.diagnostics",
    "category": "function",
    "text": "diagnostics(ss::StateSpace{T}) where T\n\nRun diagnostics and print results\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.residuals",
    "page": "Manual",
    "title": "StateSpaceModels.residuals",
    "category": "function",
    "text": "residuals(ss::StateSpace{T}) where T\n\nObtain standardized residuals from a state-space model\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.jarquebera",
    "page": "Manual",
    "title": "StateSpaceModels.jarquebera",
    "category": "function",
    "text": "jarquebera(e::Matrix{T}) where T\n\nRun Jarque-Bera normality test and return p-values\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.ljungbox",
    "page": "Manual",
    "title": "StateSpaceModels.ljungbox",
    "category": "function",
    "text": "ljungbox(e::Matrix{T}; maxlag::Int = 20) where T\n\nRun Ljung-Box independence test and return p-values\n\n\n\n\n\n"
},

{
    "location": "manual/#StateSpaceModels.homoscedast",
    "page": "Manual",
    "title": "StateSpaceModels.homoscedast",
    "category": "function",
    "text": "homoscedast(e::Matrix{T}) where T\n\nRun homoscedasticity test and return p-values\n\n\n\n\n\n"
},

{
    "location": "manual/#Diagnostics-1",
    "page": "Manual",
    "title": "Diagnostics",
    "category": "section",
    "text": "After the estimation is completed, diagnostics can be run over the residuals. The implemented diagnostics include the Jarque-Bera normality test, the Ljung-Box independence test, and a homoscedasticity test.diagnostics\nStateSpaceModels.residuals\nStateSpaceModels.jarquebera\nStateSpaceModels.ljungbox\nStateSpaceModels.homoscedast"
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": ""
},

{
    "location": "examples/#Air-Passengers-1",
    "page": "Examples",
    "title": "Air Passengers",
    "category": "section",
    "text": "Let\'s take the classical Air Passenger time series as an example. In order to avoid multiplicative effects, we use the well-known approach of taking the log of the series. The code is in the example folder.using CSV, StateSpaceModels, Plots, Statistics, Dates\n\n# Load the AirPassengers dataset\nAP = CSV.read(\"AirPassengers.csv\")\n\n# Take the log of the series\nlogAP = log.(Vector{Float64}(AP[:Passengers]))\n\n# Plot the data\np1 = plot(AP[:Date], logAP, label = \"Log-airline passengers\", legend = :topleft, color = :black)(Image: Log of Air Passengers time series)First we need to specify a state-space model. In this case, we\'ll utilize the basic structural model.# Create structural model with seasonality of 12 months\nmodel = structural(logAP, 12)Estimating the model gives us the trend and seasonal components of the time series.# Estimate a StateSpace structure\nss = statespace(model)\n\n# Analyze its decomposition in trend and seasonal\np2 = plot(AP[:Date], [ss.smoother.alpha[:, 1] ss.smoother.alpha[:, 3]], layout = (2, 1),\n            label = [\"Trend component\" \"Seasonal component\"], legend = :topleft)(Image: Trend and seasonal components for log of Air Passengers)We can also forecast this time series. In this example, we will forecast 24 months ahead.# Forecast 24 months ahead\nN = 24\npred, dist = forecast(ss, N)\n\n# Define forecasting dates\nfirstdate = AP[:Date][end] + Month(1)\nnewdates = collect(firstdate:Month(1):firstdate + Month(N - 1))\n\np3 = plot!(p1, newdates, pred, label = \"Forecast\")(Image: Forecast for log of Air Passengers)"
},

{
    "location": "examples/#Vehicle-tracking-1",
    "page": "Examples",
    "title": "Vehicle tracking",
    "category": "section",
    "text": "In order to illustrate one application that does not fall into any of the predefined models, thus requiring a user-defined model, let us consider an example from control theory. More precisely, we are going to use StateSpaceModels.jl to track a vehicle from noisy sensor data. In this case, y_t is a 2 times 1 observation vector representing the corrupted measurements of the vehicle\'s position on the two-dimensional plane in instant t. Since sensors collect the observations with the presence of additive Gaussian noise, we need to filter the observation in order to obtain a better estimate of the vehicle\'s position. The full code to run this example is in the example folder.The position and speed in each dimension compose the state of the vehicle. Let us refer to x_t^(d) as the position on the axis d and to dotx^(d)_t as the speed on the axis d in instant t. Additionally, let eta^(d)_t be the input drive force on the axis d, which acts as state noise. For a single dimension, we can describe the vehicle dynamics asbeginalign\n    x_t+1^(d) = x_t^(d) + Big( 1 - fracrho Delta_t2 Big) Delta_t dotx^(d)_t + fracDelta^2_t2 eta_t^(d) \n    dotx^(d)_t+1 = (1 - rho) dotx^(d)_t + Delta_t eta^(d)_t\nendalignwhere Delta_t is the time step and rho is a known damping effect on speed. We can cast this dynamical system as a state-space model in the following manner:beginalign \n    y_t = beginbmatrix 1  0  0  0  0  0  1  0 endbmatrix alpha_t+1 + varepsilon_t \n    alpha_t+1 = beginbmatrix 1  (1 - tfracrho Delta_t2) Delta_t  0  0  0  (1 - rho)  0  0  0  0  1  (1 - tfracrho Delta_t2)  0  0  0  (1 - rho) endbmatrix alpha_t + beginbmatrix tfracDelta^2_t2  0  Delta_t  0  0  tfracDelta^2_t2  0  Delta_t endbmatrix eta_t\nendalignwhere alpha_t = (x_t^(1) dotx^(1)_t x_t^(2) dotx^(2)_t)^top and eta_t = (eta^(1)_t eta^(2)_t)^top.We can formulate the vehicle tracking problem in StateSpaceModels.jl as:# State transition matrix\nT = kron(Matrix{Float64}(I, p, p), [1 (1 - ρ * Δ / 2) * Δ; 0 (1 - ρ * Δ)])\n# Input matrix\nR = kron(Matrix{Float64}(I, p, p), [.5 * Δ^2; Δ])\n# Output (measurement) matrix\nZ = kron(Matrix{Float64}(I, p, p), [1 0])\n# User defined model\nmodel = StateSpaceModel(y, Z, T, R)\n# Estimate vehicle speed and position\nss = statespace(model)In this example, we define the noise variances H and Q, generate the noises and simulate a random vehicle trajectory using the state-space equations:# Generate random actuators\nQ = .5 * Matrix{Float64}(I, q, q)\nη = MvNormal(Q)\n# Generate random measurement noise\nH = 2. * Matrix{Float64}(I, p, p)\nε = MvNormal(H)\n# Simulate vehicle trajectory\nα = zeros(n + 1, m)\ny = zeros(n, p)\nfor t in 1:n\n    y[t, :] = Z * α[t, :] + rand(ε)\n    α[t + 1, :] = T * α[t, :] + R * rand(η)  \nendAn illustration of the results can be seen in the following figure. It can be seen that the measurements are reasonably noisy when compared to the true position. Furthermore, the estimated positions, represented by the smoothed state, effectively estimate the true positions with small inaccuracies.(Image: Vehicle tracking)"
},

{
    "location": "reference/#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": ""
},

{
    "location": "reference/#Base.size",
    "page": "Reference",
    "title": "Base.size",
    "category": "function",
    "text": "size(model::StateSpaceModel)\n\nReturn the dimensions n, p, m and r of the StateSpaceModel\n\n\n\n\n\n"
},

{
    "location": "reference/#StateSpaceModels.ztr",
    "page": "Reference",
    "title": "StateSpaceModels.ztr",
    "category": "function",
    "text": "ztr(model::StateSpaceModel)\n\nReturn the state space model arrays Z, T and R of the StateSpaceModel\n\n\n\n\n\n"
},

{
    "location": "reference/#StateSpaceModels.check_steady_state",
    "page": "Reference",
    "title": "StateSpaceModels.check_steady_state",
    "category": "function",
    "text": "check_steady_state(P_t1::Matrix{T}, P_t::Matrix{T}, tol::T) where T\n\nReturn true if steady state was attained with respect to tolerance tol, false otherwise. The steady state is checked by the following equation maximum(abs.((P_t1 - P_t)./P_t1)) < tol.\n\n\n\n\n\n"
},

{
    "location": "reference/#StateSpaceModels.ensure_pos_sym!",
    "page": "Reference",
    "title": "StateSpaceModels.ensure_pos_sym!",
    "category": "function",
    "text": "ensure_pos_sym!(M::Matrix{T}; ϵ::T = 1e-8) where T\n\nEnsure that matrix M is positive and symmetric to avoid numerical errors when numbers are small by doing (M + M\')/2 + ϵ*I\n\n\n\n\n\n"
},

{
    "location": "reference/#StateSpaceModels.prepare_forecast",
    "page": "Reference",
    "title": "StateSpaceModels.prepare_forecast",
    "category": "function",
    "text": "prepare_forecast(ss::StateSpace{Typ}, N::Int) where Typ\n\nAdjust matrix Z for forecasting and check for dimension errors.\n\n\n\n\n\n"
},

{
    "location": "reference/#Reference-1",
    "page": "Reference",
    "title": "Reference",
    "category": "section",
    "text": "size\nztr\nStateSpaceModels.check_steady_state\nStateSpaceModels.ensure_pos_sym!\nStateSpaceModels.prepare_forecast"
},

]}
