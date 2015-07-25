require("utils")
lb = require("lb")

NX = 128
NY = 64

PARAMETERS = lb.LBGKParameters()

PARAMETERS.tau = 0.8

PY = 0

VT = 0.0 
VB = 0.0

START = 0 
END = 2048
FREQ = 128

cmap = lb.colormap()
cmap:append_color({1, 0, 0})
cmap:append_color({1, 1, 1})
cmap:append_color({0, 0, 1})

poiseuilleVelocity = function(iY, Ly, U) 
	y = iY/Ly
	return 4.*U*(y-y*y)
end

poiseuillePressure = function(iX, Lx, Ly, U, nu)
	Lx = Lx - 1
	Ly = Ly - 1
	return 8.*nu*U/(Ly*Ly)*(Lx/2.-iX);
end

density2rgb = function(rho, ux, uy)
	return cmap:map_value((rho - 2.7)/3)
end

initializer = function(x, y)
	return 2.5 
--	return 2.5 + math.random(), 0, 0
end

initialize = function(simulation, initializer, pinfo)
	--local pinfo = lb.LBPartitionInfo()
	--simulation:partition_info(pinfo)
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			local ux = poiseuilleVelocity(y, ny-1, 0.01)
			simulation:set_averages(x, y, 0,
						initializer(x + x0, y + y0), 
						ux, 
						0.0)
		end
	end

	simulation:set_equilibrium()
end

-- ##### Begin simulation code #####

simulation = lb.d2q9_BGK(1, NX, NY, 1, PY)
pinfo = lb.LBPartitionInfo()
simulation:partition_info(pinfo)
--math.randomseed(pinfo:processor_rank() + 1)

--pinfo.size[0] = 1

--print(pinfo.size[0])

simulation:set_parameters(PARAMETERS)

--[[if not PY then
	simulation:set_walls_speed(VT, VB)
end--]]

initialize(simulation, initializer, pinfo)

t0 = lb.wtime()
last_report = t0

--simulation:dump("init.h5")

for t = START, END do
	simulation:advance()

	if pinfo:processor_rank() == 0 then
		local now = lb.wtime()
		if now - last_report > 7 then
			report_progress(t0, START, END, t)
			last_report = now		
		end
	end

	if math.fmod(t, FREQ) == 0 then
		render(density2rgb, simulation, pinfo, t, ".png")
--		simulation:dump(make_filename(t))
		print(simulation:mass())
	end
end

par = lb.LBGKParameters()
simulation:get_parameters(par)

print("TAU = ", par.tau)
