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

--------------------------------------------------------------------------------------------
-- tof function (from TofDaqR)
--------------------------------------------------------------------------------------------
function M.tof(arg) 

  local toftype = arg.toftype or "LTOF"
  local drift = arg.drift or 6000
  local pulse = arg.pulse or 1000
  local massToCharge = arg.massToCharge or 100
  local x = arg.x or 0
  local v = arg.v or 0

  local amu = 1.660538921e-27  -- atomic mass unit (kg)
  local e = 1.60217657e-19  -- elementary charge (C)

  local Vdrift = -abs(drift)
  local Vpush = abs(pulse)
  local Vpull = -Vpush

  local d1 = 0.004 + 0.0035  -- pull to push distance

  local u1 = Vpush - Vpull
  local u3 = Vpull - Vdrift
  local d3, d4, d5, d6, d7
  local u5, u6, u7

  if (toftype == "HTOF") then
    d3 = 0.0065
    d4 = 0.507 + 0.511
    d5 = 0.0165
    d6 = 0.0713 - d5
  elseif (toftype == "HTOF-W") then
    d3 = 0.0065
    d4 = 0.507 + 0.511 + 2 * 0.5125
    d5 = 0.0165
    d6 = 0.0713 - d5
    d7 = 0.0085
    u7 = 1.5 * Vpush - Vdrift
  elseif (toftype == "LTOF") then
    d3 = 0.014
    d4 = 1.0465 + 1.052
    d5 = 0.0385
    d6 = 0.1593 - d5
  elseif (toftype == "CTOF") then
    d3 = 0.006
    d4 = 2 * 0.1435
    d5 = 0.017
    d6 = 0.0165
  elseif (toftype == "NTOF") then
    d3 = 0.014
    d4 = 3*0.8
    d5 = 0.0385
    d6 = 0.1593 - d5
    d7 = d5
    u7 = 1.5*Vpush - Vdrift
  else
    error("invalid toftype")
  end

  -- calculate u5 and u6
  local x0 = d1-0.001
  local k0 = x0/d1
  local p0 = k0+u3/u1
  local a, b
  if (toftype == "HTOF-W") then
    a = d1/u1*k0^(-0.5) + d3/u3*(p0^(-0.5) - k0^(-0.5)) - d4/u1/2*p0^(-1.5) + 2*d7/u7*p0^(-0.5)
    b = d1/u1/2*k0^(-1.5) + d3/u3/2*(p0^(-1.5) - k0^(-1.5)) - d4/u1*3/4*p0^(-2.5) + d7/u7*p0^(-1.5)
    u5 = (a-2*p0*b + 4*d5/u1*p0^(-1.5))/(-2*b/u1)
    u6 = (-4*d6*(p0-u5/u1)^(-0.5))/(a+4*d5/u5*(p0^(-0.5) - (p0-u5/u1)^(-0.5)))
  elseif (toftype == "NTOF") then
    a = d1/u1*k0^(-0.5) + d3/u3*(p0^(-0.5) - k0^(-0.5)) - d4/u1/2*p0^(-1.5) + 2*d7/u7*p0^(-0.5)
    b = d1/u1/2*k0^(-1.5) + d3/u3/2*(p0^(-1.5) - k0^(-1.5)) - d4/u1*3/4*p0^(-2.5) + d7/u7*p0^(-1.5)
    u5 = (a-2*p0*b + 2*d5/u1*p0^(-1.5))/(-2*b/u1)
    u6 = (-2*d6*(p0-u5/u1)^(-0.5))/(a+2*d5/u5*(p0^(-0.5) - (p0-u5/u1)^(-0.5)))
  else
    a = d1/u1*k0^(-0.5) + d3/u3*(p0^(-0.5) - k0^(-0.5)) - d4/u1/2*p0^(-1.5)
    b = d1/u1/2*k0^(-1.5) + d3/u3/2*(p0^(-1.5) - k0^(-1.5)) - d4/u1*3/4*p0^(-2.5)
    u5 = (a-2*p0*b + 2*d5/u1*p0^(-1.5))/(-2*b/u1)
    u6 = (-2*d6*(p0-u5/u1)^(-0.5))/(a+2*d5/u5*(p0^(-0.5) - (p0-u5/u1)^(-0.5)))
  end

  local xi = v * sqrt(massToCharge*amu / (2 * e*u1))
  local k = (x0 - x) / d1

  local timeOfFlight

  -- calculate time-of-flight
  if (toftype == "HTOF-W") then
    timeOfFlight = sqrt(massToCharge*amu/(2*e)) *
      (2*d1/sqrt(u1)*(sqrt(xi*xi + k) - xi) +
      2*d3/u3*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1)) +
      d4/sqrt(u1*xi*xi + k*u1 + u3) +
      2*4*d5/u5*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1 + u3 - u5)) +
      2*4*d6/u6*sqrt(u1*xi*xi + k*u1 + u3 - u5) +
      4*d7/u7*sqrt(u1*xi*xi + k*u1 + u3))
  elseif (toftype == "NTOF") then
    timeOfFlight = sqrt(massToCharge*amu/(2*e)) *
      (2*d1/sqrt(u1)*(sqrt(xi*xi + k) - xi) +
      2*d3/u3*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1)) +
      d4/sqrt(u1*xi*xi + k*u1 + u3) +
      4*d5/u5*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1 + u3 - u5)) +
      4*d6/u6*sqrt(u1*xi*xi + k*u1 + u3 - u5) +
      4*d7/u7*sqrt(u1*xi*xi + k*u1 + u3))
  else
    timeOfFlight = sqrt(massToCharge*amu/(2*e)) *
      (2*d1/sqrt(u1)*(sqrt(xi*xi + k) - xi) +
      2*d3/u3*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1)) +
      d4/sqrt(u1*xi*xi + k*u1 + u3) +
      4*d5/u5*(sqrt(u1*xi*xi + k*u1 + u3) - sqrt(u1*xi*xi + k*u1 + u3 - u5)) +
      4*d6/u6*sqrt(u1*xi*xi + k*u1 + u3 - u5))
  end

  return timeOfFlight
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
