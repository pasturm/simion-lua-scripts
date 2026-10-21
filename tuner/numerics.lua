--[[
numerics.lua: helper for tuner.lua

Householder-QR least squares, quadratic response-surface model,
bounded derivative-free optimization on [-1,1]^n.

Patrick Sturm
(c) TOFWERK, 2026
--]]

local M = {}

local function finite(x)
  return type(x) == 'number' and x == x and x ~= math.huge and x ~= -math.huge
end
M.finite = finite

local function clamp(x, lo, hi)
  if x < lo then return lo elseif x > hi then return hi else return x end
end

local function acopy(x)           -- array (shallow) copy
  local y = {}
  for i=1,#x do y[i] = x[i] end
  return y
end

-- least squares ---------------------------------------------------------- 
function M.least_squares(A, b, tolerance)
  local m = #A
  assert(m > 0 and #b == m, 'least-squares data mismatch')
  local n = #A[1]
  assert(m >= n, 'too few rows for least-squares model')
  local R, y = {}, {}
  for i=1,m do
    assert(#A[i] == n, 'ragged model matrix')
    R[i] = {}
    for j=1,n do R[i][j] = A[i][j] end
    y[i] = b[i]
  end
  local scale = 0
  for i=1,m do for j=1,n do scale = math.max(scale, math.abs(R[i][j])) end end
  local tol = (tolerance or 1e-12)*math.max(1, scale)*math.max(m, n)
  for k=1,n do
    local norm = 0
    for i=k,m do norm = norm+R[i][k]*R[i][k] end
    norm = math.sqrt(norm)
    assert(norm > tol, 'rank-deficient response-surface model')
    local sign = R[k][k] >= 0 and 1 or -1
    local v = {}
    for i=k,m do v[i] = R[i][k] end
    v[k] = v[k]+sign*norm
    local vv = 0
    for i=k,m do vv = vv+v[i]*v[i] end
    assert(vv > 0, 'QR decomposition failed')
    local beta = 2/vv
    for j=k,n do
      local dot = 0
      for i=k,m do dot = dot+v[i]*R[i][j] end
      dot = beta*dot
      for i=k,m do R[i][j] = R[i][j]-dot*v[i] end
    end
    local dot = 0
    for i=k,m do dot = dot+v[i]*y[i] end
    dot = beta*dot
    for i=k,m do y[i] = y[i]-dot*v[i] end
  end
  local x = {}
  for i=n,1,-1 do
    local s = y[i]
    for j=i+1,n do s = s-R[i][j]*x[j] end
    assert(math.abs(R[i][i])>tol, 'rank-deficient response-surface model')
    x[i] = s/R[i][i]
    assert(finite(x[i]), 'non-finite least-squares coefficient')
  end
  return x
end

-- quadratic RSM (per response) ------------------------------------------------
function M.feature_count(n) return 1 + 2*n + n*(n-1)/2 end

function M.features(x)
  local z = {1}
  for i=1,#x do z[#z+1] = x[i] end
  for i=1,#x-1 do for j=i+1,#x do z[#z+1] = x[i]*x[j] end end
  for i=1,#x do z[#z+1] = x[i]*x[i] end
  return z
end

function M.fit(coded_rows, response_values)
  assert(#coded_rows == #response_values, 'RSM data length mismatch')
  local X,y = {},{}
  for i=1,#coded_rows do
    if finite(response_values[i]) then     -- drop failed/non-finite runs
      X[#X+1] = M.features(coded_rows[i])
      y[#y+1] = response_values[i]
    end
  end
  local p = M.feature_count(#coded_rows[1])
  assert(#X >= p, 'too few valid response rows for full second-order model')
  local beta = M.least_squares(X,y)
  local model = {coefficients = beta, nfactor = #coded_rows[1], n = #X}
  function model:predict(x)
    local f = M.features(x)
    local v = 0
    for i=1,#f do v = v+self.coefficients[i]*f[i] end
    return v
  end
  return model
end

-- deterministic bounded optimizer -----------------------------------------
local primes = {2,3,5,7,11,13,17}
local function halton(index, base)
  local f,r = 1,0
  while index > 0 do
    f = f/base
    r = r+f*(index%base)
    index = math.floor(index/base)
  end
  return r
end

local function refine(fn, start, step, tol, max_iter)
  local x = acopy(start)
  local fx = fn(x)
  local it = 0
  while step > tol and it < max_iter do
    it = it+1
    local improved = false
    for j=1,#x do
      for _,sgn in ipairs({-1,1}) do
        local z = acopy(x)
        z[j] = clamp(z[j]+sgn*step, -1, 1)
        local fz = fn(z)
        if fz < fx then x,fx,improved = z,fz,true end
      end
    end
    if not improved then step = step/2 end
  end
  return x,fx
end

function M.minimize(fn, n, opt)
  opt = opt or {}
  assert(n >= 1 and n <= 7, 'optimizer supports 1..7 dimensions')
  local samples = opt.samples or 257
  local starts = opt.starts or 12
  local step = opt.step or 0.25
  local tol = opt.tol or 1e-6
  local maxit = opt.max_iter or 3000
  local cand = {}
  local center = {}
  for j=1,n do center[j] = 0 end
  cand[1] = {x = center, f = fn(center)}
  for i=1,samples do
    local x = {}
    for j=1,n do x[j] = 2*halton(i,primes[j])-1 end
    cand[#cand+1] = {x = x, f = fn(x)}
  end
  table.sort(cand, function(a,b) return a.f < b.f end)
  local count = math.min(starts, #cand)
  local bx,bf = nil,math.huge
  for i=1,count do
    local x,f = refine(fn, cand[i].x, step, tol, maxit)
    if f < bf then bx,bf = x,f end
  end
  return bx,bf
end

function M.maximize(fn, n, opt)
  local x,f = M.minimize(function(z) return -fn(z) end, n, opt)
  return x,-f
end

return M