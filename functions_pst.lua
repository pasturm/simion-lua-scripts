--[[
functions_pst.lua
Statistics and other utility functions

Patrick Sturm
(c) 2017-2020 TOFWERK
 --]]
 
local M = {}

-------------------------------------------------------------------------------------------------------
-- get kinetic energy of particle
-------------------------------------------------------------------------------------------------------
function M.get_ke()
	local speed = math.sqrt(ion_vx_mm^2 + ion_vy_mm^2 + ion_vz_mm^2)
	return speed_to_ke(speed, ion_mass)
	-- ke = 1/2*ion_mass*1.660538921e-27*(speed*1000)^2 / 1.602176565e-19
end

-------------------------------------------------------------------------------------------------------
-- Some statistics functions
-------------------------------------------------------------------------------------------------------
-- Get the mean value of a table
function M.mean(t)
	local sum = 0
	local count = 0
	for _,v in pairs(t) do
		if type(v) == 'number' and v == v then  -- strip non-numbers and nans
			sum = sum + v
			count = count + 1
		end
	end
	return  count ~= 0 and sum/count or nil
end

-- Get the standard deviation of a table
function M.sd(t)
	local m = M.mean(t)
	local vm
	local sum = 0
	local count = 0
	for _,v in pairs(t) do
		if type(v) == 'number' and v == v then  -- strip non-numbers and nans
			vm = v - m
			sum = sum + (vm * vm)
			count = count + 1
		end
	end
	return count > 1 and (math.sqrt(sum / (count-1))) or nil
end

-- Get quantiles of a table (p = 0.25: 1st quartile, p = 0.5: median, p = 0.75: 3rd quartile).
function M.quantile(t,p)
	local temp={}

	-- deep copy table so that when we sort it, the original is unchanged
	-- also weed out any non-numbers and nans
	for _,v in pairs(t) do
		if type(v) == 'number' and v == v then
			table.insert(temp, v)
		end
	end
	table.sort(temp)

	local h
	local quant = 0
	if #temp > 1 then
	    -- this is based on R's default method (type 7, see help(quantile) in R)
		h = ((#temp - 1) * p) + 1
		-- h = ((#temp + 1/3) * p) + 1/3  -- R's method 8 (probably best method, see help(quantile) in R)
		quant = temp[floor(h)] + ((h - floor(h)) * (temp[floor(h) + 1] - temp[floor(h)]))
	end
	return quant ~= 0 and quant or nil
end

-- Get the interquartile range of a table
function M.iqr(t)
	local intqr = nil
	if #t > 1 then
		intqr = M.quantile(t,0.75) - M.quantile(t,0.25)
	end
	return intqr
end

-- Get the maximum value of a table
function M.max(t)
	local lim = -math.huge
	for _,v in pairs(t) do
		if type(v) == 'number' then
			if v > lim then
				lim = v
			end
		end
	end
	return lim
end

-- Get the minimum value of a table
function M.min(t)
	local lim = math.huge
	for _,v in pairs(t) do
		if type(v) == 'number' then
			if v < lim then
				lim = v
			end
		end
	end
	return lim
end

-- Get the root mean square value of a table
function M.rms(t)
	local sum = 0
	local count = 0
	for _,v in pairs(t) do
		if type(v) == 'number' and v == v then
			sum = sum + v * v
			count = count + 1
		end
	end
	return count ~= 0 and math.sqrt(sum/count) or nil
end

-- Get the mean absolute deviation of a table
function M.mad(t)
	local sum = 0
	local count = 0
	local m = M.mean(t)
	for _,v in pairs(t) do
		if type(v) == 'number' and v == v then
			sum = sum + math.abs(v - m)
			count = count + 1
		end
	end
	return count ~= 0 and sum/count or nil
end

--------------------------------------------------------------------------------------------
-- Split function (from http://lua-users.org/wiki/SplitJoin)
--------------------------------------------------------------------------------------------
function string:split(sSeparator, nMax, bRegexp)
	assert(sSeparator ~= '')
	assert(nMax == nil or nMax >= 1)

	local aRecord = {}

	if self:len() > 0 then
		local bPlain = not bRegexp
		nMax = nMax or -1

		local nField = 1 nStart = 1
		local nFirst,nLast = self:find(sSeparator, nStart, bPlain)
		while nFirst and nMax ~= 0 do
			aRecord[nField] = self:sub(nStart, nFirst - 1)
			nField = nField + 1
			nStart = nLast + 1
			nFirst,nLast = self:find(sSeparator, nStart, bPlain)
			nMax = nMax - 1
		end
		aRecord[nField] = self:sub(nStart)
	end

	return aRecord
end

-------------------------------------------------------------------------------------------------------
-- sign of a number
-------------------------------------------------------------------------------------------------------
function M.sign(x)
	assert(type(x) == 'number', "argument is not a number")
	return x > 0 and 1 or x < 0 and -1 or 0
end

--------------------------------------------------------------------------------------------
return M
