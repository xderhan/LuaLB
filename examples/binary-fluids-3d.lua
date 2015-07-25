lb = require("lb")

--
--  parameters
--

NX = 64
NY = 64
NZ = 64

PARAMETERS = lb.LBMixParameters()

PARAMETERS.T = 0.511 -- 0.498, 0.511, 0.526
PARAMETERS.a = 9/49
PARAMETERS.b = 2/21
PARAMETERS.K = 0.01
PARAMETERS.Gr = 1.0
PARAMETERS.lamda = 1.1
PARAMETERS.rtau = 0.8
PARAMETERS.ptau = 0.8

PZ = 1

START = 0
END = 1024
FREQ = 64

--
--  initializer
--
poiseuilleVelocity = function(iZ, Lz, U) 
	z = iZ/Lz
	return 4.*U*(z-z*z)
end

initializer = function(x, y, z)
	return 2.5 + math.random(), 0, 0, 0
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
					initializer(x + x0, y + y0, z + z0), 
					poiseuilleVelocity(z, nz-1, 0.01), 
					0.0, 
					0.0,
					0.0)
			end
		end
	end

	simulation:set_equilibrium()
end

make_filename = function(t)
	local res = ""
	local digits = 1

	for n = 1, 8 do
		pp = math.pow(10, n)
		if math.floor(t / pp) == 0 then
			break
		else
			digits = digits + 1
		end
	end

	todo = 7 - digits
	for n = 0, todo do
		res = res .. "0"
	end

	return res .. t .. ".png"
end

render = function(simulation, t)
	local filename = make_filename(t)
	local pinfo = simulation:partition_info()

	simulation:dump(string.format("tmp-%d.h5", t))

--	if pinfo.processor_rank == 0 then
--		local cmd = "./lb-iso -o " .. filename

--		os.execute(cmd)
--		os.execute("rm -f tmp.h5")
--	end
end

--
--  usefull function
--

report_progress = function(start_wtime, t0, t1, t)
	local secs = lb.wtime() - start_wtime
	local steps = t - t0
	local sps = secs/steps

	print()
	io.write(string.format("<>> %d steps done in %f seconds (%f sec/step)\n",
			       steps, secs, sps))
	io.write(string.format("<>> ETA: %s\n", lb.wtime_string((t1 - t)*sps)))
end

--
--  finally ... core functionality
--

simulation = lb.d3q19_Mix(NX, NY, NZ, PZ)

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

	local M = simulation:mass()
	print(M)

	if math.fmod(t, FREQ) == 0 then
		render(simulation, t)
	end
end
