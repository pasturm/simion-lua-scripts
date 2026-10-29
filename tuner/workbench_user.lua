--[[
  SIMION Lua workbench program
  Patrick Sturm
  (c) TOFWERK, 2026
--]]

simion.workbench_program()

local TUNE = simion.import('tuner.lua')
local FILE = simion.import('fileio.lua')

adjustable _V1 = 1.374599
adjustable _V2 = 1.115826
adjustable _V3 = 0.744575
adjustable _V4 = 0
adjustable _V5 = -0.898785

adjustable _tune_enabled = 0  -- voltage tuning
adjustable _load_best = 0     -- load last optimized voltages

local resolution   -- TOF resolution
local sensitivity  -- TOF sensitivity

-- Tuner configuration.
local config = {}
config.write = false
config.stop = {
  max_repeats = 10,
  rel_tol     = 1e-3,
  abs_tol     = 0,
}
config.responses = {
  { Name = "res",  Target = 1e6, Weight = 1 },
  { Name = "sens", Target = 1,   Weight = 0 },
}
-- Between 3 and 7 variables must be enabled.
config.variables = {
  { Name = "V1", Enabled = 1, Range = 0.005, LowLimit = 0, HighLimit = 2 },
  { Name = "V2", Enabled = 1, Range = 0.005, LowLimit = 0.5, HighLimit = 1.5 },
  { Name = "V3", Enabled = 1, Range = 0.005, LowLimit = 0.5, HighLimit = 1.5 },
  { Name = "V4", Enabled = 0, Range = 0.02, LowLimit = -0.5, HighLimit = 0.5 },
  { Name = "V5", Enabled = 1, Range = 0.005, LowLimit = -1.5, HighLimit = -0.5 },
}

local function current_values()
  return {V1 = _V1, V2 = _V2, V3 = _V3, V4 = _V4, V5 = _V5}
end

local hooks = {}

function hooks.apply_variables(v)
  _V1 = v.V1
  _V2 = v.V2
  _V3 = v.V3
  _V4 = v.V4
  _V5 = v.V5
end

function hooks.run()
  run()
  return { res  = resolution, sens = sensitivity }
end

-- function hooks.prerun(run, run_once)
-- end

function segment.flym()
  if _tune_enabled == 1 then
    result = TUNE.run(config, hooks, current_values(), _load_best)
  else
    run()
  end
end

function segment.init_p_values()
  if _load_best == 0 or _tune_enabled == 1 then
    adj_elect01 = _V1
    adj_elect02 = _V2
    adj_elect03 = _V3
    adj_elect04 = _V4
    adj_elect05 = _V5
  else
    local best = assert(FILE.load_best('tuning_results'),
      'no saved best point found in tuning_results')
    adj_elect01 = best.variables.V1
    adj_elect02 = best.variables.V2
    adj_elect03 = best.variables.V3
    adj_elect04 = best.variables.V4
    adj_elect05 = best.variables.V5
  end
end

function segment.terminate_run()
  -- Response variables are configured here:
  -- resolution = ...
  -- sensitivity = ...
end