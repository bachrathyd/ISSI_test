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

function CompRand(x::Vector)
    return rand(size(x)...) .- 0.5 + 1.0im*(rand(size(x)...) .- 0.5)
end
function CompRand(x::Real)
    return rand() .- 0.5 + 1.0im*(rand() .- 0.5)
end
function CompRand(x)
    return rand() .- 0.5 + 1.0im*(rand() .- 0.5)
end

function dynamic_problem(prob,alg,maxdelay;Historyresolution=20,eigN=8,zerofixpont=true,p=[]);
   
    StateSmaplingTime=LinRange(-maxdelay,0.0,Historyresolution); 
    #first solution for random initial functions
    #TODO: the true random thing, leads to too small stap size at the initial simulation
    solset= [
        solve(
            remake(prob;u0=CompRand(u0), h=(p,t)->u0+CompRand(u0))
            ,alg) for k in 1:eigN];

    #StateSmaplingTimeSMALL=LinRange(-maxdelay,0.0,100);  
    #TODO: melyik kezdő függvény a jó?!?!?!?! - lehet sima lépcsős is?

  #  interpolationfunctions=[
  #           LinearInterpolation(StateSmaplingTime,
  #           [[cos((k)*t),sin((k)*t)] for t in StateSmaplingTime]
  #           ) for k in 1:eigN
  #           ]

  #interpolationfunctions=[
  #  LinearInterpolation(StateSmaplingTime,
  #   [[cos((k)*t)] for t in StateSmaplingTime]
  #     ) for k in 1:eigN
  #         ];

  #interpolationfunctions=[
  #  LinearInterpolation(StateSmaplingTime,
  #   [[(prob.tspan[end]+t<prob.tspan[end]*k/eigN && prob.tspan[end]+t>prob.tspan[end]*(k-1)/eigN) ? 1.0 : 0.0] for t in StateSmaplingTime]
 #      ) for k in 1:eigN
  #         ]         ;
    
  #              solset= [
  #                  solve(
  #                      remake(prob;u0=[2.0+1im,0.5+1im]*0.0 +1.0*CompRand(u0),#p=(k,prob.p[2:end]...),
  #                      h=(p,t)->[2.0+1im,0.5+1im] .*(k*1+0) .* interpolationfunctions[k](t))
  #                      ,alg,reltol=1e-4,abstol=1e-4) for k in 1:eigN];

       #solset= [
       #                     solve(
       #                         remake(prob;u0=1.0 .+ 0.0*CompRand(u0),#p=(k,prob.p[2:end]...),
       #                         h=(p,t)->[exp(-1.0im*k*t)*2.0,exp(-1.0im*k*t)])
        #                        ,alg,reltol=1e-5,abstol=1e-5) for k in 1:eigN];
#
    #solset= [
    #    solve(
     #       remake(prob;u0=CompRand(u0),
     #       h=(p,t)->LinearInterpolation(StateSmaplingTimeSMALL,[[sin(k*t/maxdelay),sin(k*t/maxdelay)] for t in StateSmaplingTimeSMALL])(t))
     #       ,alg,reltol=1e-3,abstol=1e-3) for k in 1:eigN];
    StateCombinations=diagm(0 => ones(typeof(CompRand(prob.u0[1])),eigN))## at the initaliatio there is no combination, a random inital values are used

    dynamic_problem(prob,alg,maxdelay,StateSmaplingTime,eigN,zerofixpont,solset,StateCombinations);

end

function iterate!(dp::dynamic_problem);
    #TODO: az S-ek alapján csinálni az új history függvényt, bár ennek nem kellene, hogy értelme legyen.
    dp_prev=deepcopy(dp);
    for k in 1:dp.eigN
        hloc(p,t)=histopryremap(t+T,dp_prev.StateCombinations,dp_prev.SolutionSet,k)
        dp.SolutionSet[k]= solve(
            remake(dp.DDEdynProblem;u0= hloc(p,0.0), 
            h=(p,t)-> hloc(p,t))
            ,alg,reltol=1e-3,abstol=1e-3)
    end
end

function compute_eig!(dp::dynamic_problem)
    #TODO: az előző végén lévő adatokból meg lehet határozni: a Si=Vi-1*dp.StateCombinations 
    Si=[getvalues(dp.SolutionSet[solind],t) for t in dp.StateSmaplingTime, solind in 1:dp.eigN];
    Vi=[getvalues(dp.SolutionSet[solind],(t+dp.DDEdynProblem.tspan[end])) for t in dp.StateSmaplingTime, solind in 1:dp.eigN];#initialization of the Starting Vector

    SV=Si'*Vi;
    SS=Si'*Si;
    
    Snorm=[norm(Si[:,k]) for k in 1:dp.eigN ];
    
    Hi=SS\SV .+ 0.0im;

    μs,Ai=eigen(Hi, sortby = x -> -abs(x));

    
    #dp.StateCombinations[:] .= ((Ai))[:]# ./ μs./ Snorm)[:];#length(dp.StateSmaplingTime)*  ./ Snorm   ./ μs
   dp.StateCombinations[:] .= ((Ai) / diagm(μs))[:];#length(dp.StateSmaplingTime)*  ./ Snorm   ./ μs
   # dp.StateCombinations = Ai;#./ μs)[:];#length(dp.StateSmaplingTime)*  ./ Snorm   ./ μs
    return  μs,Si,Vi,Ai
end
#TODO: check saveat=....
#In this case it saves only the values at the necessary timeponts? BUT will it be still prcise enough for the next integrations????
function spectralRadiusOfMapping(dp::dynamic_problem)
    #TODO: itt valami rendesebb iteráció kell, vagy akár a gyökök számát is autómatikusan változtatni
    for k=1:10
        compute_eig!(dp);
        iterate!(dp);
    end
    ei,si,vi,aii=compute_eig!(dp);
    spectraradiuse=maximum(abs.(ei));
    @show spectraradiuse
    return spectraradiuse,ei
end

function histopryremap(t,A::Matrix,SolutionSetV,solindex)
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


function SVi1real(dp::dynamic_problem,idx)
    Si=[(getvalues(dp.SolutionSet[solind],t))[idx] for t in dp.StateSmaplingTime, solind in 1:dp.eigN];
    Vi=[(getvalues(dp.SolutionSet[solind],(t+dp.DDEdynProblem.tspan[end])))[idx] for t in dp.StateSmaplingTime, solind in 1:dp.eigN];#initialization of the Starting Vector
 return  Si,Vi
end
