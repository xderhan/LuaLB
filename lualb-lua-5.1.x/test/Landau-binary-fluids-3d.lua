dofile("utils.lua")

--
--  parameters
--

NX = 64
NY = 64
NZ = 64

PARAMETERS = lb.LBLandauParameters()

PARAMETERS.a = -0.1
PARAMETERS.b = 0.1
PARAMETERS.K = 0.09
PARAMETERS.rtau = 1.0
PARAMETERS.ptau = 1.0

PZ = 1

START = 0
END = 2048
FREQ = 64

VB = {0,0}
VT = {0,0}

--
--  initializer
--

initializer = function(x, y, z)
	return 1+0.0001*(1-2*math.random()-1.), 0.0001*(1-2*math.random()), 0, 0, 0
end

--
--  initializes averages
--

initialize = function(simulation, initializer)
	local pinfo = simulation:partition_info()
	local nx, ny, nz = pinfo:size(0), pinfo:size(1), pinfo:size(2)
	local x0, y0, z0 = pinfo:global_origin(0),
			   pinfo:global_origin(1),
			   pinfo:global_origin(2)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			for z = 0, nz - 1 do
				simulation:set_averages(x, y, z,
					initializer(x + x0, y + y0, z + z0))
			end
		end
	end

	simulation:set_equilibrium()
end

render3 = function(simulation, t)
	local filename = make_filename(t,".h5")
	local pinfo = simulation:partition_info()

	simulation:dump(filename)

--	if pinfo.processor_rank == 0 then
--		local cmd = "./lb-iso -o " .. filename

--		os.execute(cmd)
--		os.execute("rm -f tmp.h5")
--	end

	collectgarbage() -- this should on C side
end

--
--  finally ... core functionality
--

simulation = lb.d3q19_LD(NX, NY, NZ, PZ)
pinfo = simulation:partition_info()
math.randomseed(pinfo:processor_rank() + 1)

simulation:set_parameters(PARAMETERS)

if not PZ then
	simulation:set_walls_speed(VT, VB)
end

initialize(simulation, initializer)

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
	
	print(simulation:mass(),"\n")

	if math.mod(t, FREQ) == 0 then
		render3(simulation, t)
	end
end
