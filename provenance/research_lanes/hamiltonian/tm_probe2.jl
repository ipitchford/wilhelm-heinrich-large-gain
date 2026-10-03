# Feasibility probe: rigorous (Taylor-model) enclosure of the time-T flow of the WH system from a crest box.
# This is a cost/width probe for PROGRAM.md, not a certificate of anything.
using TaylorModels, TaylorModels.ValidatedInteg, IntervalArithmetic, TaylorSeries, Printf

function wh!(dx, x, p, t)
    mu, k = p
    dx[1] = x[2]
    dx[2] = k - mu*x[2] - x[3]
    dx[3] = k*exp(x[1]) - x[3]
    nothing
end

function probe(mu, k, U0, y0, T, delta; orderQ=2, orderT=18, abstol=1e-14, maxsteps=40000, which=1)
    q0 = [U0, 0.0, y0]
    X0 = [interval(q0[1]-delta, q0[1]+delta), interval(-1e-9, 1e-9), interval(q0[3]-delta*k, q0[3]+delta*k)]
    sp = JetSpace(2*orderQ, ["a","b","c"])
    local sol
    el = @elapsed begin
        if which == 1
            sol = validated_integ(wh!, X0, 0.0, T, orderQ, orderT, abstol, sp, (mu, k); maxsteps=maxsteps, parse_eqs=false)
        else
            sol = validated_integ2(wh!, X0, 0.0, T, orderQ, orderT, abstol, sp, (mu, k); maxsteps=maxsteps, parse_eqs=false)
        end
    end
    tt = expansion_point(sol)
    n = length(sol)
    qTM = sol[n]
    domt = domain(qTM[1])
    tend = tt[end] + sup(domt)
    dtI = interval(sup(domt), sup(domt))
    nb = [interval(-1.0, 1.0) for _ in 1:3]
    q = [evaluate(evaluate(qTM[i], dtI), nb) for i in 1:3]
    w = [diam(q[i]) for i in 1:3]
    @printf("mu=%g k=%g T=%.3f delta=%.0e orderQ=%d orderT=%d integ%d: steps=%d reached t=%.3f (%.0f%%) time=%.2fs  set widths at t_end U=%.2e p=%.2e y=%.2e\n",
            mu, k, T, delta, orderQ, orderT, which, n, tend, 100*tend/T, el, w[1], w[2], w[3])
    flush(stdout)
    return sol
end

cases = [ (0.05, 0.07875, 3.0528722958389216, 1.1444732475352755, 56.66331745737882),
          (0.05, 0.1575, 4.67809378669197, 7.004299811467587, 112.13291332447197) ]

println("# warm-up / compile")
probe(1.0, 2.4, 1.3004, 4.9363, 0.5, 1e-4)
for (mu,k,U0,y0,T) in cases
    for delta in (1e-4, 1e-3, 1e-2)
        try
            probe(mu, k, U0, y0, T, delta)
        catch e
            println("FAILED mu=$mu k=$k delta=$delta: ", sprint(showerror, e)[1:min(end,200)])
        end
    end
end
