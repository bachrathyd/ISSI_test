
function histopryremapp(t,A::Matrix,SolutionSetV,solindex)
    #t [-τ,0]
    sum([SolutionSetV[i](t) * A[i,solindex] for i in 1:length(SolutionSetV)])
end


#function extract_solutionpoints(sol,time)
#    [solset[solind](t) for t in timelocations, solind in 1:eigvecnum];
#end