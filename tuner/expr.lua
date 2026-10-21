--[[
expr.lua: helper for tuner.lua

-- Arithmetic expression evaluator for factor/control transformations.
-- Supported: numbers, variables, parentheses, unary +/- and + - * / ^.

Patrick Sturm
(c) TOFWERK, 2026
--]]

local M = {}

local function tokenize(s)
  local t, i, n = {}, 1, #s
  while i <= n do
    local c = s:sub(i,i)
    if c:match('%s') then
      i = i + 1
    elseif c:match('[%a_]') then
      local j = i + 1
      while j <= n and s:sub(j,j):match('[%w_]') do j = j + 1 end
      t[#t+1] = {kind = 'id', value = s:sub(i,j-1)}
      i = j
    elseif c:match('[%d.]') then
      local rest = s:sub(i)
      local num = rest:match('^%d+%.?%d*[eE][%+%-]?%d+')
        or rest:match('^%.%d+[eE][%+%-]?%d+')
        or rest:match('^%d+%.?%d*')
        or rest:match('^%.%d+')
      assert(num and tonumber(num), 'invalid number near: '..rest)
      t[#t+1] = {kind = 'num', value = tonumber(num)}
      i = i + #num
    elseif c:match('[%+%-%*/%^%(%)]') then
      t[#t+1] = {kind = c, value = c}
      i = i + 1
    else
      error('unsupported character in expression: '..c)
    end
  end
  t[#t+1] = {kind = 'eof'}
  return t
end

function M.eval(expression, env)
  if type(expression) == 'function' then return expression(env) end
  assert(type(expression) == 'string', 'transformation must be a string or function')
  local tok, p = tokenize(expression), 1
  local parse_expression
  local function peek() return tok[p] end
  local function take(k)
    local x = tok[p]
    assert(x.kind == k, "expected '"..k.."', got '"..x.kind.."'")
    p = p + 1
    return x
  end
  local function primary()
    local x = peek()
    if x.kind == 'num' then p = p+1
      return x.value
    end
    if x.kind == 'id' then
      p = p+1
      local v = env[x.value]
      assert(type(v) == 'number', 'unknown/non-numeric variable: '..x.value)
      return v
    end
    if x.kind == '(' then
      p = p+1
      local v = parse_expression()
      take(')')
      return v
    end
    error('expected number, variable, or parenthesis')
  end
  local function unary()
    if peek().kind == '+' then 
      p = p+1
      return unary() 
    end
    if peek().kind == '-' then p = p+1
      return -unary() 
    end
    return primary()
  end
  local function power()
    local a = unary()
    if peek().kind == '^' then p = p+1
      return a ^ power()
    end
    return a
  end
  local function product()
    local a = power()
    while peek().kind == '*' or peek().kind == '/' do
      local op = peek().kind
      p = p+1
      local b = power()
      if op == '*' then a = a*b else a=a/b end
    end
    return a
  end
  parse_expression = function()
    local a = product()
    while peek().kind == '+' or peek().kind == '-' do
      local op = peek().kind
      p = p+1
      local b = product()
      if op == '+' then a = a+b else a = a-b end
    end
    return a
  end
  local value = parse_expression()
  take('eof')
  assert(value == value and value ~= math.huge and value ~= -math.huge,
    'transformation returned a non-finite value')
  return value
end

return M
