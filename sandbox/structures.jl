struct dynamic_problem
    DDEdynProblem::DDEProblem #e.g. =DDEProblem(....)
    alg::MethodOfSteps # e.g.: alg = MethodOfSteps(Tsit5())
    maxdelay::Float64
    StateSmaplingTime::Vector
    eigN::Integer # number of eigen vectors
    zerofixpont::Bool
    SolutionSet::Vector{ODESolution}
    StateCombinations::Matrix
    #EigenVectors::Matrix
    #EigenValues::Matrix
end

function dynamic_problem(prob,alg,maxdelay;Historyresolution=20,eigN=8,zerofixpont=true,reltol=1e-3,abstol=1e-3);
        
    StateSmaplingTime=-maxdelay  :maxdelay/(Historyresolution-1) :0;
    #first solution for random initial functions
    #TODO: the true random thing, leads to too small stap size at the initial simulation
    solset= [
        solve(
            remake(prob;u0=CompRand(u0), h=(p,t)->CompRand(u0))
            ,alg,reltol=1e-3,abstol=1e-3) for k in 1:eigN];

    StateCombinations=diagm(0 => ones(typeof(CompRand(prob.u0[1])),eigN))## at the initaliatio there is no combination, a random inital values are used

    dynamic_problem(prob,alg,maxdelay,StateSmaplingTime,eigN,zerofixpont,solset,StateCombinations);

end

function iterate!(dp::dynamic_problem);
    dp_prev=deepcopy(dp);
    for k in 1:dp.eigN
        hloc(p,t)=histopryremapp(t+T,dp_prev.StateCombinations,dp_prev.SolutionSet,k)
        dp.SolutionSet[k]= solve(
            remake(prob;u0= hloc(p,0), 
            h=(p,t)-> hloc(p,t))
            ,alg,reltol=1e-3,abstol=1e-3)
    end
end

function CompRand(x::Vector)
    return rand(size(x)...) .- 0.5 + 1.0im*(rand(size(x)...) .- 0.5)
end
function CompRand(x::Real)
    return rand() .- 0.5 + 1.0im*(rand() .- 0.5)
end

function compute_eig!(dp::dynamic_problem)
    #TODO: az előző végén lévő adatokból meg lehet határozni: a Si=Vi-1*dp.StateCombinations 
    Si=[getvalues(dp.SolutionSet[solind],t) for t in dp.StateSmaplingTime, solind in 1:dp.eigN];
    Vi=[getvalues(dp.SolutionSet[solind],(t+dp.DDEdynProblem.tspan[end])) for t in dp.StateSmaplingTime, solind in 1:dp.eigN];#initialization of the Starting Vector

    SV=Si'*Vi;
    SS=Si'*Si;
    
    Snorm=[norm(Si[:,k]) for k in 1:dp.eigN ];
    
    Hi=SS\SV;

    μs,Ai=eigen(Hi);
    μs,Ai=sorteigen(μs,Ai);
    dp.StateCombinations[:] .= (Ai ./ μs)[:];#length(dp.StateSmaplingTime)*  ./ Snorm   ./ μs
    return  μs
end


function SVi1real(dp::dynamic_problem,idx)
    Si=[(getvalues(dp.SolutionSet[solind],t))[idx] for t in dp.StateSmaplingTime, solind in 1:dp.eigN];
    Vi=[(getvalues(dp.SolutionSet[solind],(t+dp.DDEdynProblem.tspan[end])))[idx] for t in dp.StateSmaplingTime, solind in 1:dp.eigN];#initialization of the Starting Vector
 return  Si,Vi
end

function histopryremapp(t,A::Matrix,SolutionSetV,solindex)
    #sum([SolutionSetV[i](t) * A[i,solindex] for i in 1:length(SolutionSetV)])
    sum([getvalues(SolutionSetV[i],t) * A[i,solindex]  for i in 1:length(SolutionSetV)])
end




function getvalues(sol::ODESolution,t::Real)
    if t<0
        sol.prob.h(sol.prob.p,t);
    elseif t==0
        sol.prob.u0;
    else
        sol(t)
    end
end

function getvalues(sol::ODESolution,t::Vector)
    [getvalues(sol,ti)[1] for ti in t]
end

function sorteigen(evals::Vector{T},evecs::Matrix{T}) where {T}
    p = sortperm(abs.(evals))
    evals[p[end:-1:1]], evecs[:, p[end:-1:1]]
end