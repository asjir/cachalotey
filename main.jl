using WAV, DSP, Statistics, Plots, Unitful
plotly()
audio, fs = wavread("terminal.wav")

wavplay(audio, fs)

audio[:, 1] |> histogram

audiolen = length(audio[:, 1]) ./ fs
plot(eachindex(audio[:, 1]) ./ fs, audio[:, 1])

# these values are arbitrary - 15-20 simply gives enough points for calculated standard deviations to be stable
window = 15
overlap = 5

"""compute standard deviations of signal over windows
 note: if I took maximum over windows then the values wouldn't be independent from the size of the window - https://en.wikipedia.org/wiki/Order_statistic"""
adjstds(audio; window, overlap) =
    map(sqrt ∘ (*), std.(arraysplit(audio, window + overlap, overlap)), std.(arraysplit(audio, window, 0)))

"""if we assume there's no pause less than say 0.01s, then any shorter interruption is part of same click"""
anyininterval(a; wps, minimalinterval) =
    map(any ∘ (x -> Bool.(x)), arraysplit(a, 1 + Int(ceil(minimalinterval * wps)), Int(ceil(minimalinterval * wps))))

plotsec(a; window) = plot(eachindex(a)u"s" ./ fs * window, a)

# looking at this plot can give idea for threshold
plotsec(adjstds(audio[:, 1]; window, overlap); window)
threshold = 0.03

# and this for minimal interval between the clicks
plotsec(adjstds(audio[:, 1]; window, overlap) .> threshold; window)
minimalinterval = 0.01

function findclicktimes(audio; window, overlap, threshold, minimalinterval)
    wps = fs / window
    isclick = adjstds(audio; window, overlap) .> threshold
    isclick2 = anyininterval(isclick; wps, minimalinterval)
    startsclick = map(x -> iszero(first(x)) & isone(last(x)), arraysplit(isclick2, 2, 1))
    clickstarttimes = (eachindex(startsclick)./wps)[startsclick]
    endsclick = map(x -> isone(first(x)) & iszero(last(x)), arraysplit(isclick2, 2, 1))
    clickendtimes = (eachindex(endsclick)./wps)[endsclick]
    (clickstarttimes .+ clickendtimes) ./ 2
end

ct_nofilter = findclicktimes(audio[:, 1]; window, overlap, threshold, minimalinterval)
ct_nofilter[2] .- ct_nofilter[1]

# https://www.imeko.org/publications/tc19-Metrosea-2020/IMEKO-TC19-MetroSea-2020-19.pdf
# gives 8000+Hz for click
# first 0.25s is noise, so let's check it out
pg = periodogram(audio[1:Int(fs * 0.25), 1]; fs)
# let's filter out 99% of power
pg.freq[findfirst(cumsum(pg.power) .> 0.99 * sum(pg.power))]
# 4200Hz way below 8kHz


hpass = digitalfilter(Highpass(4000; fs), Butterworth(4))
filtered = filt(hpass, audio[:, 1])
plotsec(adjstds(filtered; window, overlap); window)
# after filtering it's expected we need a new threshold
threshold = 0.008

ct_filter = findclicktimes(filtered; window, overlap, threshold, minimalinterval)

ct_filter
ct_nofilter

# interestingly applying filter classified the last sound aat 3.59 as not a coda

# instead of eyeballing thresholds we can plot results for different values

let idxs = 0.01 * 10 .^ (-1:0.1:1), nclicksfound(threshold) = length(findclicktimes(filtered; window, overlap, threshold, minimalinterval))
    plot(idxs, nclicksfound.(idxs); xscale=:log10)
end

let idxs = 0.01 * 10 .^ (-1:0.1:1), nclicksfound(threshold) = length(findclicktimes(audio[:, 1]; window, overlap, threshold, minimalinterval))
    plot(idxs, nclicksfound.(idxs); xscale=:log10)
end

# this can be used to quickly inspect which clicks are detected correctly and which not
first_differing = findfirst((!isapprox).(findclicktimes(filtered; window, overlap, threshold=1e-2, minimalinterval), findclicktimes(filtered; window, overlap, threshold=5e-3, minimalinterval)[1:end-1]; atol=minimalinterval / 2))

findclicktimes(filtered; window, overlap, threshold=1e-2, minimalinterval)[first_differing], findclicktimes(filtered; window, overlap, threshold=5e-3, minimalinterval)[first_differing]
# 0.685s, 0.663s

# but same procedure for minimalinterval makes sense only for the highpass-filtered input
let idxs = 0.01 * 10 .^ (-1:0.1:1), nclicksfound(minimalinterval) = length(findclicktimes(filtered; window, overlap, threshold, minimalinterval))
    plot(idxs, nclicksfound.(idxs); xscale=:log10)
end

let idxs = 0.01 * 10 .^ (-1:0.1:1), nclicksfound(minimalinterval) = length(findclicktimes(audio[:, 1]; window, overlap, threshold, minimalinterval))
    plot(idxs, nclicksfound.(idxs); xscale=:log10)
end

clickintervals(clicktimes) = clicktimes[2:end] .- clicktimes[1:end-1]

histogram(clickintervals(ct_filter) .* 750u"m"; bins=100)

ci_filter = clickintervals(ct_filter)

Plots.scatter(ci_filter[1:end-1]u"s", ci_filter[2:end]u"s")