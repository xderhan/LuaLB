dofile("utils.lua")

NX = 128
NY = 128

PARAMETERS = lb.LBLandauParameters()

PARAMETERS.a = -0.1
PARAMETERS.b = 0.1
PARAMETERS.K = 0.09
PARAMETERS.rtau = 1.0
PARAMETERS.ptau = 1.0

PY = 1

VT = 0.05 
VB = -0.05 

START = 0 
END = 2048 -- 1024, 2048, 4096, 8192, 16384, 32768
FREQ = 128

cmap = lb.colormap()
cmap:append_color({1, 0, 0})
cmap:append_color({1, 1, 1})
cmap:append_color({0, 0, 1})

density2rgb = function(rho, ux, uy)
	return cmap:map_value(rho)
end

ordering2rgb = function(phi, ux, uy)	
	return cmap:map_value( (1+phi)/2 )
end

initializer = function(x, y)
	return 1+0.0001*(1-2*math.random()-1.), 0.0001*(1-2*math.random()), 0, 0
end

initialize = function(simulation, initializer)
	local pinfo = simulation:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			simulation:set_averages(x, y,
						initializer(x + x0, y + y0))
		end
	end

	simulation:set_equilibrium()
end

-- ##### Begin simulation code #####

simulation = lb.d2q9_LD(NX, NY, PY)
pinfo = simulation:partition_info()
math.randomseed(pinfo:processor_rank() + 1)

simulation:set_parameters(PARAMETERS)

if not PY then
	simulation:set_walls_speed(VT, VB)
end

initialize(simulation, initializer)
--simulation:dump("bin_init_2d.h5")

t0 = lb.wtime()
last_report = t0

for t = START, END do
	simulation:advance()

	if pinfo:processor_rank() == 0 then
		local now = lb.wtime()
		if now - last_report > 7 then
			report_progress(t0, START, END, t)
			last_report = now
		end
	end

	stats = simulation:stats()
	print("ENERGY\t= ", stats.kin_energy)
	print("MASS\t= ", simulation:mass())

	if math.mod(t, FREQ) == 0 then
		render(ordering2rgb, simulation, t)
--		simulation:dump(make_filename(t, ".h5"))
	end
end

--simulation:dump("final.h5")

simulation:destroy()
