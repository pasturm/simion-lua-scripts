--[[
fileio.lua: helper for tuner.lua

Result, best-point and resume handling.

Patrick Sturm
(c) TOFWERK, 2026
--]]

local M = {}

local LFS = require "lfs"  -- Load LuaFileSystem library.

local function serialize(v,indent)
  indent = indent or ''
  if type(v) == 'number' or type(v) == 'boolean' then return tostring(v)
  elseif type(v) == 'string' then return string.format('%q',v)
  elseif type(v) == 'table' then
    local keys = {}
    for k in pairs(v) do keys[#keys+1] = k end
    table.sort(keys,function(a,b) return tostring(a) < tostring(b) end)
    local out = {'{'}
    for _,k in ipairs(keys) do
      local key = type(k) == 'number' and ('['..k..']') or ('['..string.format('%q',k)..']')
      out[#out+1] = '\n'..indent..'  '..key..'='..serialize(v[k],indent..'  ')..','
    end
    out[#out+1] = '\n'..indent..'}'
    return table.concat(out)
  elseif v == nil then return 'nil' end
  error('cannot serialize '..type(v))
end

local function save_lua(path,value)
  local tmp = path..'.tmp'
  local f = assert(io.open(tmp, 'w'))
  f:write('return ', serialize(value), '\n')
  f:close()
  os.remove(path)
  assert(os.rename(tmp, path))
end

local function load_lua(path)
  local f = io.open(path, 'r')
  if not f then return nil end
  f:close()
  local chunk,err = loadfile(path)
  assert(chunk, err)
  return chunk()
end

local function quote(v)
  if v == nil then return '' end
  local s = tostring(v)
  if s:find('[,"\r\n]') then return '"'..s:gsub('"','""')..'"' end
  return s
end

function M.load_best(dir)
  local root = load_lua(dir.."\\"..'bestpoint.lua')
  if root then return root end
  local newest,best = nil,nil
  local ok,iter,state,var = pcall(LFS.dir,dir)
  if not ok then return nil end
  for name in iter,state,var do
    if name ~= '.' and name ~= '..' then
      local p = dir.."\\"..name.."\\"..'bestpoint.lua'
      local a = LFS.attributes(p)
      if a and a.mode == 'file' and (not newest or a.modification > newest) then
        local value = load_lua(p)
        if value then newest,best = a.modification,value end
      end
    end
  end
  return best
end

function M.save_best(dir, best)
  LFS.mkdir(dir)
  save_lua(dir.."\\"..'bestpoint.lua', best)
  local f = assert(io.open(dir.."\\"..'bestpoint.txt', 'w'))
  f:write('Name|Value\n')
  local keys = {}
  for k in pairs(best.controls or {}) do 
    keys[#keys+1] = k 
  end
  table.sort(keys)
  for _,k in ipairs(keys) do f:write(k,'|',tostring(best.controls[k]),'\n') end
  f:close()
end

function M.save_repeat(dir, repeat_no, records,best)
  LFS.mkdir(dir)
  local sub = dir.."\\"..os.date('%Y-%m-%d-%Hh%Mm%Ss')..string.format('_%02d',repeat_no)
  LFS.mkdir(sub)
  local f = assert(io.open(sub.."\\"..'results.csv','w'))
  local factor_keys,control_keys,response_keys = {},{},{}
  if records[1] then
    for k in pairs(records[1].factors or {}) do factor_keys[#factor_keys+1]=k end
    for k in pairs(records[1].controls or {}) do control_keys[#control_keys+1]=k end
    for k in pairs(records[1].raw or {}) do response_keys[#response_keys+1]=k end
  end
  table.sort(factor_keys)
  table.sort(control_keys)
  table.sort(response_keys)
  local h = {'kind','repeat','run','measured_desirability'}
  for _,k in ipairs(factor_keys) do h[#h+1] = 'factor_'..k end
  for _,k in ipairs(control_keys) do h[#h+1] = 'control_'..k end
  for _,k in ipairs(response_keys) do h[#h+1] = 'response_'..k end
  f:write(table.concat(h,','),'\n')
  for _,r in ipairs(records) do
    local row = {r.kind, r.repeat_index, r.run_index, r.measured_desirability}
    for _,k in ipairs(factor_keys) do row[#row+1] = r.factors[k] end
    for _,k in ipairs(control_keys) do row[#row+1] = r.controls[k] end
    for _,k in ipairs(response_keys) do row[#row+1] = r.raw[k] end
    for i=1,#row do row[i] = quote(row[i]) end
    f:write(table.concat(row,','),'\n')
  end
  f:close()
  save_lua(sub.."\\"..'bestpoint.lua',best)
  local tf = assert(io.open(sub.."\\"..'bestpoint.txt','w'))
  tf:write('Name|Value\n')
  local keys = {}
  for k in pairs(best.controls or {}) do 
    keys[#keys+1] = k 
  end
  table.sort(keys)
  for _,k in ipairs(keys) do 
    tf:write(k,'|', tostring(best.controls[k]), '\n') 
  end
  tf:close()
end

return M
