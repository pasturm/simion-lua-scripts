--[[
tuner.lua: DOE/RSM voltage tuner for SIMION.

Box-Behnken design -> 
independent full quadratic RSM per response -> 
deterministic optimization of combined desirability ->
verification run -> 
recenter on best measured point -> 
repeat.

Patrick Sturm
(c) TOFWERK, 2026
--]]

local EXPR = simion.import('expr.lua')
local NUM  = simion.import('numerics.lua')
local FILE = simion.import('fileio.lua')

local output_dir = 'tuning_results'
local M = {}

local function copy(t)  -- deep copy (config/records)
  local r = {}
  for k,v in pairs(t or {}) do r[k] = type(v) == 'table' and copy(v) or v end
  return r
end

-- Box-Behnken design --------------------------------------------
local bb_blocks = {
  [3] = {{1,2},{1,3},{2,3}},
  [4] = {{1,2},{3,4},{1,4},{2,3},{1,3},{2,4}},
  [5] = {{1,2},{1,3},{3,4},{4,5},{2,5},{1,4},{1,5},{2,3},{2,4},{3,5}},
  [6] = {{1,2,4},{2,3,5},{3,4,6},{1,4,5},{2,5,6},{1,3,6}},
  [7] = {{4,5,6},{1,6,7},{2,5,7},{1,2,4},{3,4,7},{1,3,5},{2,3,6}},
}
local function box_behnken(k)
  assert(bb_blocks[k], 'Box-Behnken design supports 3..7 enabled factors')
  local out, m = {}, #bb_blocks[k][1]
  local n = 2^m
  for _,b in ipairs(bb_blocks[k]) do
    for r=0,n-1 do
      local x = {}
      for j=1,k do x[j] = 0 end
      for j=1,m do x[b[j]] = (math.floor(r/2^(j-1))%2==0) and -1 or 1 end
      out[#out+1] = x
    end
  end
  local z = {}
  for j=1,k do z[j] = 0 end
  out[#out+1] = z  -- single center point
  return out
end

-- desirability ----------------------------------------
local function desirability(y, target, pmin, pmax)  -- lower is better
  local P = y < target and pmin or pmax
  local denom = target - P
  if denom == 0 then return 0 end
  local X = (target - y)/denom
  return 0.2*X + 0.8*X*X - 1
end
local function desirability_overall(values, responses, pmins, pmaxs)
  local sw,total = 0,0
  for i,r in ipairs(responses) do
    local w = r.Weight or 1
    total = total + w*desirability(values[i], r.Target, pmins[i], pmaxs[i])
    sw = sw + w
  end
  assert(sw ~= 0, 'sum of response weights must not be zero')
  return total/sw
end

-- factor <-> control transforms ----------------------------------------
local function enabled_factors(cfg)
  local e = {}
  for _,f in ipairs(cfg.factors) do
    if f.Enabled == true or f.Enabled == 1 then e[#e+1] = f end
  end
  assert(#e >= 3 and #e <= 7, 'between 3 and 7 factors must be enabled')
  return e
end
local function controls_env(cfg)
  local env = {}
  for _,c in ipairs(cfg.controls) do env[c.Name] = c.StartValue end
  return env
end
local function derive_centers(cfg, controls)
  local centers = {}
  for _,f in ipairs(cfg.factors) do
    local v = EXPR.eval(f.Transformation, controls)
    if f.Enabled == true or f.Enabled == 1 then
      local lo,hi = f.LowLimit+f.Range, f.HighLimit-f.Range
      assert(lo <= hi, 'factor range does not fit limits: '..f.Name)
      if v < lo then v = lo elseif v > hi then v = hi end
    end
    centers[f.Name] = v
  end
  return centers
end
local function decode(coded, enabled, centers)
  local physical = copy(centers)
  for i,f in ipairs(enabled) do physical[f.Name] = centers[f.Name]+f.Range*coded[i] end
  return physical
end
local function make_controls(cfg, factors)
  local c = {}
  for _,def in ipairs(cfg.controls) do c[def.Name] = EXPR.eval(def.Transformation, factors) end
  return c
end
local function set_start_values(cfg, controls)
  for _,c in ipairs(cfg.controls) do
    if controls[c.Name] ~= nil then c.StartValue = controls[c.Name] end
  end
end

-- run execution ----------------------------------------------------
local function execute_run(hooks, run)
  assert(type(hooks.apply_controls) == 'function','hooks.apply_controls required')
  assert(type(hooks.run) == 'function','hooks.run required')
  hooks.apply_controls(run.controls, run)
  local in_prerun = false
  local function run_once()
    assert(in_prerun, 'run_once is valid only inside hooks.prerun')
    return hooks.run(run)
  end
  if hooks.prerun then
    assert(type(hooks.prerun) == 'function', 'hooks.prerun must be a function')
    in_prerun = true
    local ok,err = pcall(hooks.prerun, run, run_once)
    in_prerun = false
    if not ok then error(err,0) end
  end
  local raw = assert(hooks.run(run), 'final run did not return a response table')
  assert(type(raw) == 'table', 'final run must return a response table')
  return copy(raw)
end
local function response_vector(raw, responses)
  local y = {}
  for i,r in ipairs(responses) do y[i] = raw[r.Name] end
  return y
end
local function run_record(kind, rep, index, coded, factors, controls, raw, responses)
  return {kind = kind, repeat_index = rep, run_index = index, coded = copy(coded),
    factors = copy(factors), controls = copy(controls), raw = copy(raw),
    responses = response_vector(raw, responses)}
end

-- stopping on individual measured responses ---------------------------
local function stop_config(cfg)
  local s = cfg.stop or {}
  return {
    max_repeats = s.max_repeats or 10,
    rel_tol = s.rel_tol or 1e-3,
    abs_tol = s.abs_tol or 0,
    tol = s.tol or {},  -- per-response {abs=,rel=} by Name
  }
end
local function response_improved(r, prev, cur, sc)
  if prev == nil or cur == nil then return false end
  local T = r.Target
  local gain = math.abs(prev - T) - math.abs(cur - T)  -- >0 => closer to target
  local o = sc.tol[r.Name] or {}
  local abs_tol = o.abs or sc.abs_tol
  local rel_tol = o.rel or sc.rel_tol
  local threshold = math.max(abs_tol, rel_tol*math.max(math.abs(T), math.abs(prev)))
  return gain > threshold
end
local function made_progress(responses, prev_raw, cur_raw, sc)
  if not prev_raw then return true end
  for _,r in ipairs(responses) do
    if response_improved(r, prev_raw[r.Name], cur_raw[r.Name], sc) then
      return true  -- any single response improved -> progress
    end
  end
  return false
end

-- main loop ----------------------------------------------------------
function M.run(cfg, hooks)
  assert(type(cfg) == 'table' and type(hooks) == 'table','config and hooks required')
  local enabled = enabled_factors(cfg)
  local design = box_behnken(#enabled)
  local sc = stop_config(cfg)
  local control_values = controls_env(cfg)
  if cfg.resume then
    local prior = FILE.load_best(output_dir)
    if prior and prior.controls then
      control_values = copy(prior.controls)
      set_start_values(cfg, control_values)
    end
  end

  local last_result, prev_best_raw, stall = nil, nil, 0

  for rep=1,sc.max_repeats do
    print("Repeat "..rep.."/"..sc.max_repeats..", experiments 1 to "..#design.." running...")
    local centers = derive_centers(cfg, control_values)
    local records, coded_rows = {}, {}
    local values = {}
    for i=1,#cfg.responses do values[i] = {} end

    -- (1) Box-Behnken experiments
    for i,coded in ipairs(design) do
      local factors  = decode(coded, enabled, centers)
      local controls = make_controls(cfg, factors)
      local run = {kind = 'design', repeat_index = rep, run_index = i,
        coded = copy(coded), factors = factors, controls = controls}
      local raw = execute_run(hooks, run)
      local rec = run_record('design', rep, i, coded, factors, controls, raw, cfg.responses)
      records[#records+1] = rec
      coded_rows[#coded_rows+1] = copy(coded)
      for r=1,#cfg.responses do values[r][#coded_rows] = rec.responses[r] end
    end

    -- (2) independent full quadratic model per response
    local models = {}
    for r=1,#cfg.responses do models[r] = NUM.fit(coded_rows, values[r]) end

    -- (3) per-response predicted range for desirability scaling
    local pmins, pmaxs = {}, {}
    for r,model in ipairs(models) do
      local objective = function(x) return model:predict(x) end
      local _,mn = NUM.minimize(objective, #enabled)
      local _,mx = NUM.maximize(objective, #enabled)
      pmins[r], pmaxs[r] = mn, mx
    end

    -- (4) deterministic optimization of combined desirability
    local optimum, predicted_cost = NUM.minimize(function(x)
      local vals = {}
      for i,m in ipairs(models) do vals[i] = m:predict(x) end
      return desirability_overall(vals, cfg.responses, pmins, pmaxs)
    end, #enabled)

    -- (5) one verification run at the predicted optimum
    local vfactors = decode(optimum, enabled, centers)
    local vcontrols = make_controls(cfg, vfactors)
    local vrun = {kind = 'verification', repeat_index = rep, run_index = #design+1,
      coded = copy(optimum), factors = vfactors, controls = vcontrols}
    local vraw = execute_run(hooks, vrun)
    local verification = run_record('verification', rep, #design+1, optimum,
      vfactors, vcontrols, vraw, cfg.responses)
    verification.predicted = {}
    for i,m in ipairs(models) do verification.predicted[i] = m:predict(optimum) end
    verification.predicted_desirability = predicted_cost
    records[#records+1] = verification

    -- (6) select best measured point (design + verification); recenter
    local best
    for _,rec in ipairs(records) do
      rec.measured_desirability = desirability_overall(rec.responses, cfg.responses, pmins, pmaxs)
      if not best or rec.measured_desirability < best.measured_desirability then
        best = rec
      end
    end
    control_values = copy(best.controls)
    set_start_values(cfg, control_values)
    local best_state = {repeat_index = rep, kind = best.kind, controls = copy(best.controls),
      factors = copy(best.factors), responses = copy(best.raw),
      measured_desirability = best.measured_desirability, stall = stall}
    FILE.save_best(output_dir, best_state)
    if cfg.write then FILE.save_repeat(output_dir, rep, records, best_state) end
    last_result = {best = best_state, verification = verification, records = records,
      models = models, Pred_min = pmins, Pred_max = pmaxs, centers = centers}
    local str = ""
    for name,value in pairs(best.controls) do str = str..name.."="..string.format('%.6f',value).."|" end
    for name,value in pairs(best.raw)      do str = str..name.."="..string.format('%.0f',value).."|" end
    print(string.sub(str,1,-2))

    -- (7) stopping test on individual measured responses (scale-invariant)
    local progressed = made_progress(cfg.responses, prev_best_raw, best.raw, sc)
    if prev_best_raw and not progressed then stall = stall+1 else stall = 0 end
    prev_best_raw = copy(best.raw)
    local stop_reason
    if stall >= 1 then
      stop_reason = string.format("no meaningful improvement in individual responses")
    elseif rep >= sc.max_repeats then
      stop_reason = "reached max_repeats ("..sc.max_repeats..")"
    end
    if stop_reason then
      last_result.stop_reason = stop_reason
      -- print("Stopping after repeat "..rep..": "..stop_reason)
      break
    end
  end
  return last_result
end

return M
