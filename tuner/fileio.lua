--[[
fileio.lua: helper for tuner.lua

Result, best-point and resume handling.

Patrick Sturm
(c) TOFWERK, 2026
--]]
 
local M = {}
local LFS = require "lfs"
 
local function serialize(value, indent)
  indent = indent or ''
  if type(value) == 'number' or type(value) == 'boolean' then
    return tostring(value)
  elseif type(value) == 'string' then
    return string.format('%q', value)
  elseif type(value) == 'table' then
    local keys = {}
    for key in pairs(value) do
      keys[#keys + 1] = key
    end
    table.sort(keys, function(a, b) return tostring(a) < tostring(b) end)
    local output = {'{'}
    for _, key in ipairs(keys) do
      local serialized_key
      if type(key) == 'number' then
        serialized_key = '[' .. key .. ']'
      else
        serialized_key = '[' .. string.format('%q', key) .. ']'
      end
      output[#output + 1] = '\n' .. indent .. '  ' .. serialized_key .. '=' ..
        serialize(value[key], indent .. '  ') .. ','
    end
    output[#output + 1] = '\n' .. indent .. '}'
    return table.concat(output)
  elseif value == nil then
    return 'nil'
  end
  error('cannot serialize ' .. type(value))
end
 
local function save_lua(path, value)
  local temporary = path .. '.tmp'
  local file = assert(io.open(temporary, 'w'))
  file:write('return ', serialize(value), '\n')
  file:close()
  os.remove(path)
  assert(os.rename(temporary, path))
end
 
local function load_lua(path)
  local file = io.open(path, 'r')
  if not file then
    return nil
  end
  file:close()
  local chunk, err = loadfile(path)
  assert(chunk, err)
  return chunk()
end
 
local function quote(value)
  if value == nil then
    return ''
  end
  local text = tostring(value)
  if text:find('[,"\r\n]') then
    return '"' .. text:gsub('"', '""') .. '"'
  end
  return text
end
 
local function save_best_text(path, best)
  local file = assert(io.open(path, 'w'))
  file:write('Name|Value\n')
  local keys = {}
  for key in pairs(best.variables or {}) do
    keys[#keys + 1] = key
  end
  table.sort(keys)
  for _, key in ipairs(keys) do
    file:write(key, '|', tostring(best.variables[key]), '\n')
  end
  file:close()
end
 
function M.load_best(dir)
  local root = load_lua(dir .. '\\' .. 'bestpoint.lua')
  if root then
    return root
  end
  local newest
  local best
  local ok, iterator, state, variable = pcall(LFS.dir, dir)
  if not ok then
    return nil
  end
  for name in iterator, state, variable do
    if name ~= '.' and name ~= '..' then
      local path = dir .. '\\' .. name .. '\\' .. 'bestpoint.lua'
      local attributes = LFS.attributes(path)
      if attributes and attributes.mode == 'file'
          and (not newest or attributes.modification > newest) then
        local value = load_lua(path)
        if value then
          newest = attributes.modification
          best = value
        end
      end
    end
  end
  return best
end
 
function M.save_best(dir, best)
  LFS.mkdir(dir)
  save_lua(dir .. '\\' .. 'bestpoint.lua', best)
  save_best_text(dir .. '\\' .. 'bestpoint.txt', best)
end
 
function M.save_repeat(dir, repeat_no, records, best)
  LFS.mkdir(dir)
  local subdir = dir .. '\\' .. os.date('%Y-%m-%d-%Hh%Mm%Ss') ..
    string.format('_%02d', repeat_no)
  LFS.mkdir(subdir)
  local file = assert(io.open(subdir .. '\\' .. 'results.csv', 'w'))
  local variable_keys = {}
  local response_keys = {}
  if records[1] then
    for key in pairs(records[1].variables or {}) do
      variable_keys[#variable_keys + 1] = key
    end
    for key in pairs(records[1].raw or {}) do
      response_keys[#response_keys + 1] = key
    end
  end
  table.sort(variable_keys)
  table.sort(response_keys)
  local header = {'kind', 'repeat', 'run', 'measured_desirability'}
  for _, key in ipairs(variable_keys) do
    header[#header + 1] = key
  end
  for _, key in ipairs(response_keys) do
    header[#header + 1] = 'response_' .. key
  end
  file:write(table.concat(header, ','), '\n')
  for _, record in ipairs(records) do
    local row = {
      record.kind,
      record.repeat_index,
      record.run_index,
      record.measured_desirability,
    }
    for _, key in ipairs(variable_keys) do
      row[#row + 1] = record.variables[key]
    end
    for _, key in ipairs(response_keys) do
      row[#row + 1] = record.raw[key]
    end
    for i = 1, #row do
      row[i] = quote(row[i])
    end
    file:write(table.concat(row, ','), '\n')
  end
  file:close()
  save_lua(subdir .. '\\' .. 'bestpoint.lua', best)
  save_best_text(subdir .. '\\' .. 'bestpoint.txt', best)
end

return M
