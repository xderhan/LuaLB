-- require("./utils.lua")

NX = 128
NY = 64

PARAMETERS = lb.LBGKParameters()

tau = 0.8

PARAMETERS:set(tau)

PY = 0

VT = 0.0 
VB = 0.0

START = 0 
END = 2048
FREQ = 128

cmap = lb.colormap()
cmap:append_color(1, 0, 0)
cmap:append_color(1, 1, 1)
cmap:append_color(0, 0, 1)

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

initialize = function(simulation, initializer)
	local pinfo = simulation:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			for c = 0, 1 do
				local ux = poiseuilleVelocity(y, ny-1, 0.01)
				simulation:set_averages(x, y, c,
						initializer(x + x0, y + y0), 
						ux, 
						0.0)
			end
		end
	end

	simulation:set_equilibrium()
end

-- ##### Begin simulation code #####
my_collide = function(simulation)
	local pinfo = simulation:partition_info()
	local nx, ny = pinfo:size(0), pinfo:size(1)
	local x0, y0 = pinfo:global_origin(0), pinfo:global_origin(1)

	for x = 0, nx - 1 do
		for y = 0, ny - 1 do
			for c = 0, 1 do
				r,ux,uy = simulation:get_averages(x, y, c)
			
				local rux = r*ux
				local ruy = r*uy
	
				local ruxx = rux*ux
				local ruyy = ruy*uy
	
				local ruxy = 2.*rux*uy	
				local rusq = ruxx+ruyy
	
				-- v0 - v3 
				simulation:collide_node(0,x,y,c, (1./9.)*(r+3.*rux+4.5*ruxx-1.5*rusq))
				simulation:collide_node(1,x,y,c, (1./9.)*(r+3.*ruy+4.5*ruyy-1.5*rusq))
				simulation:collide_node(2,x,y,c, (1./9.)*(r-3.*rux+4.5*ruxx-1.5*rusq))
				simulation:collide_node(3,x,y,c, (1./9.)*(r-3.*ruy+4.5*ruyy-1.5*rusq))

				-- v4 - v7 
				simulation:collide_node(4,x,y,c, (1./36.)*(r+3.*( rux+ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq))
				simulation:collide_node(5,x,y,c, (1./36.)*(r+3.*(-rux+ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq))
				simulation:collide_node(6,x,y,c, (1./36.)*(r+3.*(-rux-ruy)+4.5*(ruxx+ruyy+ruxy)-1.5*rusq))
				simulation:collide_node(7,x,y,c, (1./36.)*(r+3.*( rux-ruy)+4.5*(ruxx+ruyy-ruxy)-1.5*rusq))
			
				-- v8
				simulation:collide_node(8,x,y,c, (4./9.)*(r-1.5*rusq));
			end
		end
	end
end

simulation = lb.d2q9_BGK(2, NX, NY, 1, PY)
pinfo = simulation:partition_info()
--math.randomseed(pinfo:processor_rank() + 1)

simulation:set_parameters(PARAMETERS)

--[[if not PY then
	simulation:set_walls_speed(VT, VB)
end--]]

initialize(simulation, initializer)

t0 = lb.wtime()
last_report = t0

for t = START, END do
	--simulation:advance()
	simulation:stream()
	simulation:average()
	my_collide(simulation)
	
	if pinfo:processor_rank() == 0 then
		local now = lb.wtime()
		if now - last_report > 7 then
			report_progress(t0, START, END, t)
			last_report = now		
		end
	end

	if math.mod(t, FREQ) == 0 then
--		render(density2rgb, simulation, t)
		simulation:dump(make_filename(t,".h5"))
		print(simulation:mass())
	end
end

simulation:destroy()
