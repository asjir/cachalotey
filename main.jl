using WAV, DSP, Statistics, Plots, Unitful
plotly()
audio, fs = wavread("terminal.wav")

wavplay(audio, fs)

# plot #1
plot(eachindex(audio[:, 1]) ./ fs, audio[:, 1], title="Audio wave plot")

# here I pick values that are arbitrary - 15-20 simply gives enough points for calculated standard deviations to be stable
window = 15
overlap = 5

"""compute standard deviations of signal over windows
 note: if I took maximum over windows then the values wouldn't be independent from the size of the window - https://en.wikipedia.org/wiki/Order_statistic"""
adjstds(audio; window, overlap) =
    map(sqrt ∘ (*), std.(arraysplit(audio, window + overlap, overlap)), std.(arraysplit(audio, window, 0)))

"""if we assume there's no pause less than say 0.01s, then any shorter interruption is part of same click"""
anyininterval(a; wps, minimalinterval) =
    map(any ∘ (x -> Bool.(x)), arraysplit(a, 1 + Int(ceil(minimalinterval * wps)), Int(ceil(minimalinterval * wps))))

plotsec(a; window, kwargs...) = plot(eachindex(a)u"s" ./ fs * window, a; kwargs...)

# looking at this plot can give idea for threshold, see plot #2
plotsec(adjstds(audio[:, 1]; window, overlap); window, title="Window standard deviation")
threshold = 0.03

# and this for minimal interval between the clicks, see plot #3
plotsec(adjstds(audio[:, 1]; window, overlap) .> threshold; window, title="Threshold met binary signal")
minimalinterval = 0.01

function findclicktimes(audio; window, overlap, threshold, minimalinterval)
    wps = fs / window
    isclick = adjstds(audio; window, overlap) .> threshold
    isclick2 = anyininterval(isclick; wps, minimalinterval)
    startsclick = map(x -> iszero(first(x)) & isone(last(x)), arraysplit(isclick2, 2, 1))
    clickstarttimes = (eachindex(startsclick)./wps)[startsclick]
    endsclick = map(x -> isone(first(x)) & iszero(last(x)), arraysplit(isclick2, 2, 1))
    clickendtimes = (eachindex(endsclick)./wps)[endsclick]
    length(clickstarttimes) == length(clickendtimes) ?
    (clickstarttimes .+ clickendtimes) ./ 2 :
    typeof(clickstarttimes)[]
end

ct_nofilter = findclicktimes(audio[:, 1]; window, overlap, threshold, minimalinterval)
# 86 clicks, the last one at 3.59541

# Now I shall try to apply a high pass filter, on google I've found
# https://www.imeko.org/publications/tc19-Metrosea-2020/IMEKO-TC19-MetroSea-2020-19.pdf
# which gives 8000+Hz for clicks
# first 0.25s is noise, so let's check out what frequencies it has:
pg = periodogram(audio[1:Int(fs * 0.25), 1]; fs)
# let's filter out 99% of power
pg.freq[findfirst(cumsum(pg.power) .> 0.99 * sum(pg.power))]
# 4200Hz is way below 8kHz, so should be fine


hpass = digitalfilter(Highpass(4000; fs), Butterworth(4))
filtered = filt(hpass, audio[:, 1])
# plot #4
plot(eachindex(filtered) ./ fs, filtered, title="Audio wave plot after highpass filter")

# see plot #5.
plotsec(adjstds(filtered; window, overlap); window, title="Window standard deviation with highpass filter")
# after filtering it's expected we need a new threshold
threshold = 0.008

ct_filter = findclicktimes(filtered; window, overlap, threshold, minimalinterval)
# 85 clicks, the last one at 3.555, i.e. with filter, the 3.59541 one is not detected as click


# instead of eyeballing thresholds we can plot results for different values, this is plot #6
let idxs = 0.01 * 10 .^ (-1:0.1:1), nclicksfound(threshold) = length(findclicktimes(audio[:, 1]; window, overlap, threshold, minimalinterval))
    plot(idxs, nclicksfound.(idxs); xscale=:log10, xlabel="threshold", ylabel="#clicks detected", title="without highpass filter")
end

# plot #7
let idxs = 0.01 * 10 .^ (-1:0.1:1), nclicksfound(threshold) = length(findclicktimes(filtered; window, overlap, threshold, minimalinterval))
    plot(idxs, nclicksfound.(idxs); xscale=:log10, xlabel="threshold", ylabel="#clicks detected", title="with highpass filter")
end

# this can be used to quickly inspect when additional click gets detected, whether it was correct
first_differing = findfirst(map((x, y) -> !isapprox(x, y; atol=minimalinterval), findclicktimes(filtered; window, overlap, threshold=1e-2, minimalinterval), findclicktimes(filtered; window, overlap, threshold=5e-3, minimalinterval)))
first_differing_left, first_differing_right = minmax(findclicktimes(filtered; window, overlap, threshold=1e-2, minimalinterval)[first_differing], findclicktimes(filtered; window, overlap, threshold=5e-3, minimalinterval)[first_differing])

# in this case the additional click was the little hump, totally not a valid click, see screenshot #1
plotsec(adjstds(filtered; window, overlap); window, xlim=(first_differing_left - minimalinterval, first_differing_right + minimalinterval))

# similarly for the interval, although it should differ from recording to recording as much as threshold, because threshold depends on signal to noise ratio.
let idxs = 0.01 * 10 .^ (-1:0.1:1), nclicksfound(minimalinterval) = length(findclicktimes(filtered; window, overlap, threshold, minimalinterval))
    plot(idxs, nclicksfound.(idxs); xscale=:log10)
end

let idxs = 0.01 * 10 .^ (-1:0.1:1), nclicksfound(minimalinterval) = length(findclicktimes(audio[:, 1]; window, overlap, threshold=0.03, minimalinterval))
    plot(idxs, nclicksfound.(idxs); xscale=:log10)
end

clickintervals(clicktimes) = clicktimes[2:end] .- clicktimes[1:end-1]

# quick histogram of inter-click intervals where x axis is what distance would produce echo of this interval in water, plot #8
ci_filter = clickintervals(ct_filter)
histogram(ci_filter .* 750u"m"; bins=100, title="Distribution of intervals in meters of echo in water")

# finally a scatterplot of consecutive intervals, plot #8
covmat = cov([ci_filter[1:end-1];; ci_filter[2:end]])
Plots.scatter(ci_filter[1:end-1]u"s", ci_filter[2:end]u"s", aspect_ratio=:equal, xlim=(0, 0.14), ylim=(0, 0.14), title="consecutive intervals\nR2 score of $(covmat[2]^2 / covmat[1] / covmat[4])")
using StatsPlots
covellipse!([mean(ci_filter), mean(ci_filter)], covmat)

