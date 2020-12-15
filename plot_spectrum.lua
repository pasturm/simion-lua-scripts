--[[
 Script to plot the spectrum (FFT) of a time series.

 Patrick Sturm
 (c) 2020 TOFWERK
 --]]


-- load statistics utility functions (used for STAT.gaussian_rand())
local STAT = require "simionx.Statistics"

-- load user functions.
local FUN = simion.import("functions_pst.lua")

-- load gnuplot library code.
local PLOT = simion.import 'gnuplotlib_pst.lua'

-- local FFT = require "luafft.lua"
local FFT = simion.import("luafft.lua")

-- read data
local tof1 = {}
local file1 = "tof1.txt"
for line in io.lines(file1) do 
	table.insert(tof1, tonumber(line))
end
local tof2 = {}
local file2 = "tof2.txt"
for line in io.lines(file2) do
	table.insert(tof2, tonumber(line))
end

local tof_max = math.max(FUN.max(tof1), FUN.max(tof2))
local tof_min = math.min(FUN.min(tof1), FUN.min(tof2))

local tob_max = 1000000
local tob_min = 1
local samplingperiod = 5000  -- us

-- histogram
local hist1 = STAT.make_histogram{data=tof1, min=ceil(tof_max-tob_max+tob_min), max=floor(tof_min+tob_max-tob_min), binsize=samplingperiod, normalize = false}
local hist2 = STAT.make_histogram{data=tof2, min=ceil(tof_max-tob_max+tob_min), max=floor(tof_min+tob_max-tob_min), binsize=samplingperiod, normalize = false}

local plotdata1 = {xlabel='tof (us)', ylabel='frequency', chart_type = "scatter_lines", plot_index = 1, gnuplot_commands = "set yrange [0:]", multiplot_rows = 3, multiplot_columns = 1}

local plotdata3 = {xlabel='tof (us)', ylabel='frequency', chart_type = "scatter_lines", plot_index = 3, gnuplot_commands = "set yrange [*:]"}
for i=1,#hist1.midpoints do
	table.insert(plotdata3, {hist1.midpoints[i], hist1.frequencies[i]-hist2.frequencies[i]})
end
plot2 = PLOT.plot(plotdata3)

-- create signal
local signal = {}
for i=1,#hist1.midpoints do
  signal[i] = complex.new(hist1.frequencies[i]-hist2.frequencies[i], 0)
end

-- zero pad to next power of 2
local n = #hist1.midpoints
local nextpow2 = math.ceil(math.log(math.abs(n))/math.log(2))
local n2 = 2^nextpow2
for i=#hist1.midpoints+1,n2 do
  signal[i] = complex.new(0, 0)
end

-- subtract the mean
-- d = d - mean(d)

-- fast fourier transformation
local spec = FFT.fft(signal, false)

-- plot spectrum
local plotdata4 = {xlabel='frequency (Hz)', ylabel='amplitude', chart_type = "scatter_lines", plot_index = 4, gnuplot_commands = "set yrange [*:]"}
for i=1,n2/2 do
	table.insert(plotdata4, {(i-1)*(1/samplingperiod/n2)*1e6, spec[i]:abs()})
end
plot3 = PLOT.plot(plotdata4)
