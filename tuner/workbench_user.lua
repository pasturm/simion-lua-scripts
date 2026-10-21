--[[
SIMION Lua workbench program

Patrick Sturm
(c) TOFWERK, 2026  
 --]]

simion.workbench_program() 

local TUNE = simion.import('tuner.lua')

adjustable _V1 = 1.374599
adjustable _V2 = 1.115826
adjustable _V3 = 0.744575
adjustable _V4 = 0
adjustable _V5 = -0.898785
adjustable _tune_enabled = 0  -- voltage tuning

local resolution -- TOF resolution
local sensitivity -- TOF sensitivity

-- tuner configuration
local config = {}
config.write = false
config.resume = true
config.stop = {
  max_repeats = 10,
  rel_tol     = 1e-3,
  abs_tol     = 0,
}
config.responses = {
  { Name = "res",   Target = 1e6, Weight = 1 },
  { Name = "sens",  Target = 1,   Weight = 0 },
}
config.factors = {  -- between 3 and 7 factors need to be enabled
  { Name = "v_1",    Enabled = 1, Transformation = "V1",   Range = 0.005, LowLimit = 0,    HighLimit = 2 },
  { Name = "v_2",    Enabled = 1, Transformation = "V2",   Range = 0.005, LowLimit = 0.5,  HighLimit = 1.5 },
  { Name = "v_3",    Enabled = 1, Transformation = "V3",   Range = 0.005, LowLimit = 0.5,  HighLimit = 1.5 },
  { Name = "v_4",    Enabled = 0, Transformation = "V4",   Range = 0.02,  LowLimit = -0.5, HighLimit = 0.5 },
  { Name = "v_5",    Enabled = 1, Transformation = "V5",   Range = 0.005, LowLimit = -1.5, HighLimit = -0.5 }
}
config.controls = {
  { Name = "V1", Transformation = "v_1",  StartValue = 1.374599 },
  { Name = "V2", Transformation = "v_2",  StartValue = 1.115826 },
  { Name = "V3", Transformation = "v_3",  StartValue = 0.744575 },
  { Name = "V4", Transformation = "v_4",  StartValue = 0 },
  { Name = "V5", Transformation = "v_5",  StartValue = -0.898785 }
}

local hooks={}
function hooks.apply_controls(v)
  _V1 = v.V1
  _V2 = v.V2
  _V3 = v.V3
  _V4 = v.V4
  _V5 = v.V5
end
function hooks.run()
  run() 
  return {res = resolution, sens = sensitivity}
end
-- function hooks.prerun(run,run_once)
-- end

function segment.flym()
	if _tune_enabled == 1 then
 		result = TUNE.run(config, hooks)
	else
		run()
	end
end

function segment.init_p_values()
  adj_elect01 = _V1
  adj_elect02 = _V2
  adj_elect03 = _V3
  adj_elect04 = _V4
  adj_elect05 = _V5
end

function segment.terminate_run()
	-- response variables are configured here:
	-- resolution = ...
	-- sensitivity = ...
end
