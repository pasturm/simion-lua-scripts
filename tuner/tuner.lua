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
 
local NUM  = simion.import('numerics.lua')
local FILE = simion.import('fileio.lua')
 
local output_dir = 'tuning_results'
local M = {}
 
local function copy(t)
  local r = {}
  for k, v in pairs(t or {}) do
    r[k] = type(v) == 'table' and copy(v) or v
  end
  return r
end
 
local function is_enabled(v)
  return v.Enabled == true or v.Enabled == 1
end
 
-- Box-Behnken design -----------------------------------------------------
local bb_blocks = {
  [3] = {{1,2}, {1,3}, {2,3}},
  [4] = {{1,2}, {3,4}, {1,4}, {2,3}, {1,3}, {2,4}},
  [5] = {{1,2}, {1,3}, {3,4}, {4,5}, {2,5}, {1,4}, {1,5}, {2,3}, {2,4}, {3,5}},
  [6] = {{1,2,4}, {2,3,5}, {3,4,6}, {1,4,5}, {2,5,6}, {1,3,6}},
  [7] = {{4,5,6}, {1,6,7}, {2,5,7}, {1,2,4}, {3,4,7}, {1,3,5}, {2,3,6}}
}
 
local function box_behnken(k)
  assert(bb_blocks[k], 'Box-Behnken design supports 3..7 enabled variables')
  local out = {}
  local m = #bb_blocks[k][1]
  local n = 2^m
  for _, block in ipairs(bb_blocks[k]) do
    for r = 0, n - 1 do
      local x = {}
      for j = 1, k do x[j] = 0 end
      for j = 1, m do
        x[block[j]] = (math.floor(r / 2^(j - 1)) % 2 == 0) and -1 or 1
      end
      out[#out + 1] = x
    end
  end
  local center = {}
  for j = 1, k do center[j] = 0 end
  out[#out + 1] = center
  return out
end
 
-- Desirability ----------------------------------------------------------
local function desirability(y, target, pmin, pmax)
  local P = y < target and pmin or pmax
  local denom = target - P
  if denom == 0 then return 0 end
  local X = (target - y) / denom
  return 0.2 * X + 0.8 * X * X - 1
end
 
local function desirability_overall(values, responses, pmins, pmaxs)
  local weight_sum = 0
  local total = 0
  for i, response in ipairs(responses) do
    local weight = response.Weight or 1
    total = total + weight * desirability(values[i], response.Target, pmins[i], pmaxs[i])
    weight_sum = weight_sum + weight
  end
  assert(weight_sum ~= 0, 'sum of response weights must not be zero')
  return total / weight_sum
end
 
-- Variables -------------------------------------------------------------
local function enabled_variables(cfg)
  assert(type(cfg.variables) == 'table', 'config.variables must be a table')
  local enabled = {}
  local names = {}
  for _, variable in ipairs(cfg.variables) do
    assert(type(variable.Name) == 'string', 
      'each variable requires a Name')
    assert(not names[variable.Name], 
      'duplicate variable name: '..variable.Name)
    assert(type(variable.Range) == 'number', 
      'variable Range must be numeric: '..variable.Name)
    assert(type(variable.LowLimit) == 'number', 
      'variable LowLimit must be numeric: '..variable.Name)
    assert(type(variable.HighLimit) == 'number', 
      'variable HighLimit must be numeric: '..variable.Name)
    assert(variable.Range >= 0, 
      'variable Range must not be negative: '..variable.Name)
    assert(variable.LowLimit <= variable.HighLimit, 
      'invalid limits for variable: '..variable.Name)
    names[variable.Name] = true
    if is_enabled(variable) then enabled[#enabled + 1] = variable end
  end
  assert(#enabled >= 3 and #enabled <= 7, 
    'between 3 and 7 variables must be enabled')
  return enabled
end
 
local function check_variable_values(cfg, values)
  assert(type(values) == 'table', 
    'current variable values must be supplied')
  for _, variable in ipairs(cfg.variables) do
    assert(type(values[variable.Name]) == 'number',
      'missing or non-numeric variable value: '..variable.Name)
  end
end
 
local function derive_centers(cfg, values)
  local centers = {}
  for _, variable in ipairs(cfg.variables) do
    local value = values[variable.Name]
    if is_enabled(variable) then
      local low_center = variable.LowLimit + variable.Range
      local high_center = variable.HighLimit - variable.Range
      assert(low_center <= high_center,
        'variable range does not fit limits: ' .. variable.Name)
      if value < low_center then
        value = low_center
      elseif value > high_center then
        value = high_center
      end
    end
    centers[variable.Name] = value
  end
  return centers
end
 
local function decode(coded, enabled, centers)
  local variables = copy(centers)
  for i, variable in ipairs(enabled) do
    variables[variable.Name] =
      centers[variable.Name] + variable.Range * coded[i]
  end
  return variables
end
 
-- Run execution ---------------------------------------------------------
local function execute_run(hooks, run)
  assert(type(hooks.apply_variables) == 'function',
    'hooks.apply_variables required')
  assert(type(hooks.run) == 'function',
    'hooks.run required')
  hooks.apply_variables(run.variables, run)
  local in_prerun = false
  local function run_once()
    assert(in_prerun,
      'run_once is valid only inside hooks.prerun')
    return hooks.run(run)
  end
  if hooks.prerun then
    assert(type(hooks.prerun) == 'function',
      'hooks.prerun must be a function')
    in_prerun = true
    local ok, err = pcall(hooks.prerun, run, run_once)
    in_prerun = false
    if not ok then error(err, 0) end
  end
  local raw = assert(hooks.run(run),
    'final run did not return a response table')
  assert(type(raw) == 'table',
    'final run must return a response table')
  return copy(raw)
end
 
local function response_vector(raw, responses)
  local values = {}
  for i, response in ipairs(responses) do
    values[i] = raw[response.Name]
  end
  return values
end
 
local function run_record(kind, repeat_index, run_index, coded,
  variables, raw, responses)
  return {
    kind = kind,
    repeat_index = repeat_index,
    run_index = run_index,
    coded = copy(coded),
    variables = copy(variables),
    raw = copy(raw),
    responses = response_vector(raw, responses),
  }
end
 
-- Stopping on individual measured responses ----------------------------
local function stop_config(cfg)
  local stop = cfg.stop or {}
  return {
    max_repeats = stop.max_repeats or 10,
    rel_tol = stop.rel_tol or 1e-3,
    abs_tol = stop.abs_tol or 0,
    tol = stop.tol or {},
  }
end
 
local function response_improved(response, previous, current, stop)
  if previous == nil or current == nil then return false end
  local target = response.Target
  local gain = math.abs(previous - target) - math.abs(current - target)
  local override = stop.tol[response.Name] or {}
  local abs_tol = override.abs or stop.abs_tol
  local rel_tol = override.rel or stop.rel_tol
  local threshold = math.max(abs_tol,
    rel_tol * math.max(math.abs(target), math.abs(previous)))
  return gain > threshold
end
 
local function made_progress(responses, previous_raw, current_raw, stop)
  if not previous_raw then return true end
  for _, response in ipairs(responses) do
    if response_improved(response,
      previous_raw[response.Name], current_raw[response.Name], stop) then
      return true
    end
  end
 
  return false
end
 
-- Main loop -------------------------------------------------------------
function M.run(cfg, hooks, current_values, load_best)
  assert(type(cfg) == 'table' and type(hooks) == 'table',
    'config and hooks required')
 
  local enabled = enabled_variables(cfg)
  local design = box_behnken(#enabled)
  local stop = stop_config(cfg)
  local variable_values
 
  if load_best == true or load_best == 1 then
    local prior = assert(FILE.load_best(output_dir),
      'no saved best point found in ' .. output_dir)
    assert(type(prior.variables) == 'table',
      'saved best point has no variables table')
    variable_values = copy(prior.variables)
  else
    variable_values = copy(current_values)
  end
 
  check_variable_values(cfg, variable_values)
 
  local last_result
  local previous_best_raw
  local stall = 0
 
  for repeat_index = 1, stop.max_repeats do
    print('Repeat ' .. repeat_index .. '/' .. stop.max_repeats ..
      ', experiments 1 to ' .. #design .. ' running...')
 
    local centers = derive_centers(cfg, variable_values)
    local records = {}
    local coded_rows = {}
    local values = {}
 
    for i = 1, #cfg.responses do values[i] = {} end
 
    -- (1) Box-Behnken experiments
    for run_index, coded in ipairs(design) do
      local variables = decode(coded, enabled, centers)
      local run = {
        kind = 'design',
        repeat_index = repeat_index,
        run_index = run_index,
        coded = copy(coded),
        variables = variables,
      }
 
      local raw = execute_run(hooks, run)
      local record = run_record('design', repeat_index, run_index,
        coded, variables, raw, cfg.responses)
 
      records[#records + 1] = record
      coded_rows[#coded_rows + 1] = copy(coded)
 
      for response_index = 1, #cfg.responses do
        values[response_index][#coded_rows] = record.responses[response_index]
      end
    end
 
    -- (2) Independent full quadratic model per response
    local models = {}
    for response_index = 1, #cfg.responses do
      models[response_index] = NUM.fit(coded_rows, values[response_index])
    end
 
    -- (3) Predicted range per response for desirability scaling
    local pmins = {}
    local pmaxs = {}
 
    for response_index, model in ipairs(models) do
      local objective = function(x) return model:predict(x) end
      local _, minimum = NUM.minimize(objective, #enabled)
      local _, maximum = NUM.maximize(objective, #enabled)
      pmins[response_index] = minimum
      pmaxs[response_index] = maximum
    end
 
    -- (4) Optimize combined desirability
    local optimum, predicted_cost = NUM.minimize(function(x)
      local predicted = {}
      for i, model in ipairs(models) do
        predicted[i] = model:predict(x)
      end
      return desirability_overall(predicted, cfg.responses, pmins, pmaxs)
    end, #enabled)
 
    -- (5) Verification run at the predicted optimum
    local verification_variables = decode(optimum, enabled, centers)
    local verification_run = {
      kind = 'verification',
      repeat_index = repeat_index,
      run_index = #design + 1,
      coded = copy(optimum),
      variables = verification_variables,
    }
 
    local verification_raw = execute_run(hooks, verification_run)
    local verification = run_record('verification', repeat_index,
      #design + 1, optimum, verification_variables,
      verification_raw, cfg.responses)
 
    verification.predicted = {}
    for i, model in ipairs(models) do
      verification.predicted[i] = model:predict(optimum)
    end
    verification.predicted_desirability = predicted_cost
    records[#records + 1] = verification
 
    -- (6) Select best measured point and recenter
    local best
 
    for _, record in ipairs(records) do
      record.measured_desirability = desirability_overall(
        record.responses, cfg.responses, pmins, pmaxs)
 
      if not best or
         record.measured_desirability < best.measured_desirability then
        best = record
      end
    end
 
    variable_values = copy(best.variables)
 
    local output = ''
    for name, value in pairs(best.variables) do
      output = output .. name .. '=' .. string.format('%.6f', value) .. '|'
    end
    for name, value in pairs(best.raw) do
      output = output .. name .. '=' .. string.format('%.0f', value) .. '|'
    end
    print(string.sub(output, 1, -2))
 
    -- (7) Stop on individual measured-response progress
    local progressed = made_progress(cfg.responses, previous_best_raw, best.raw, stop)
 
    if previous_best_raw and not progressed then
      stall = stall + 1
    else
      stall = 0
    end
 
    previous_best_raw = copy(best.raw)

    local best_state = {
      repeat_index = repeat_index,
      kind = best.kind,
      variables = copy(best.variables),
      responses = copy(best.raw),
      measured_desirability = best.measured_desirability,
      stall = stall,
    }
 
    FILE.save_best(output_dir, best_state)
    if cfg.write then
      FILE.save_repeat(output_dir, repeat_index, records, best_state)
    end
 
    last_result = {
      best = best_state,
      verification = verification,
      records = records,
      models = models,
      Pred_min = pmins,
      Pred_max = pmaxs,
      centers = centers,
    }
 
    local stop_reason
    if stall >= 1 then
      stop_reason = 'no meaningful improvement in individual responses'
    elseif repeat_index >= stop.max_repeats then
      stop_reason = 'reached max_repeats (' .. stop.max_repeats .. ')'
    end
 
    if stop_reason then
      last_result.stop_reason = stop_reason
      break
    end
  end
 
  return last_result
end
 
return M
  