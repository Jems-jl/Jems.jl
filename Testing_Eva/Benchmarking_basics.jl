
# run couple of basic test scripts to fully understand everything the benchmarking do_remesh


# QUICK START TUTORIAL

using BenchmarkTools

# The `setup` expression is run once per sample, and is not included in the
# timing results. Note that each sample can require multiple evaluations
# benchmark kernel evaluations. See the BenchmarkTools manual for details.

@benchmark sort(data) setup=(data=rand(10))

##

# for quick sanity checks --> @btime

@btime sin(x) setup=(x=rand())

##

@bprofile sin(x) setup=(x=rand())

##

# If the expression you want to benchmark depends on external variables,
# you should use $ to "interpolate" them into the benchmark expression to
# avoid the problems of benchmarking with globals. Essentially, any
# interpolated variable $x or expression $(...) is "pre-computed"
# before benchmarking begins:

A = rand(3,3);

##

@btime inv($A)              # interpolate the global variable A
@btime inv($(rand(3,3)))    # rand(3,3) call occurs before benchmarking
@btime inv(rand(3,3))       # rand(3,3) call is included in benchmarking

# result --> last one takes the longest

##

# Sometimes, interpolating variables into very simple expressions
# can give the compiler more information than you intended, causing
# it to "cheat" the benchmark by hoisting the calculation out of the
# benchmark code

a = 1
b = 2

##

@btime $a * $b
@btime 1 * 2
@btime $(Ref(a))[] * $(Ref(b))[]

##


# MANUAL

# To quickly benchmark a Julia expression, use @benchmark:

@benchmark sin(1)

##

# The @benchmark macro is essentially shorthand for defining a benchmark,
# auto-tuning the benchmark's configuration parameters, and running the
# benchmark. These three steps can be done explicitly using @benchmarkable,
# tune! and run:

b = @benchmarkable sin(1);

tune!(b)        # find the right evals/sample and number of samples
                # to take for this benchmark

run(b)

## 

# @btime prints the minimum time and memory allocation before returning
# the value of the expression

@btime sin(1)

##

# @belapsed returns the minimum time in seconds.

@belapsed sin(1)

##

t = @benchmark (rand(10, 10))

##

dump(t)

##

# The minimum is a robust estimator for the location parameter of the time
# distribution, and should not be considered an outlier

minimum(t)

# The median, as a robust measure of central tendency, should be relatively
# unaffected by outliers

maximum(t)

# The mean, as a non-robust measure of central tendency, will usually be
# positively skewed by outliers

median(t)

# The maximum should be considered a primarily noise-driven outlier, and
# can change drastically between benchmark trials.

mean(t)
std(t)

##

using BenchmarkPlots, StatsPlots
b = @benchmarkable inv(rand(10,10))
t = run(b)

StatsPlots.plot(t)

##

# VISUALIZING BENCHMARK results

io = IOContext(stdout, :histmin=>0.5, :histmax=>8, :logbins=>true)
# IOContext(Base.TTY(RawFD(13) open, 0 bytes waiting))

##

b = @benchmark x^3   setup=(x = rand()); show(io, MIME("text/plain"), b)

##

b = @benchmark x^3.0 setup=(x = rand()); show(io, MIME("text/plain"), b)
BenchmarkTools.save("test.json", b)

##

d =  @benchmark x^4.0 setup=(x = rand()); show(io, MIME("text/plain"), b)
BenchmarkTools.save("test0.json", d)


##

c = BenchmarkTools.load("test.json")

##

BenchmarkTools.judge(minimum(d), minimum(b))